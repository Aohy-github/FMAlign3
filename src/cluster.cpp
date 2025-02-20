//
// Created by aohy on 24-12-25.
//

#include "cluster.h"
#include <stack>
#include "fcgr.h"
namespace Cluster {
    std::vector<int> Clu::cluster_for_lens(int longest_id) {
        int limit = All_seqs[longest_id].lens * 0.9;
        std::vector<int> result;
        All_seqs[longest_id].lens = -1;
        result.push_back(longest_id);
        for (int i = 0; i < All_seqs.size(); i++) {
            if (All_seqs[i].lens >= limit) {
                result.push_back(i);
                All_seqs[i].lens = -1;
            }
        }
        return result;
    }

    int Clu::selectLongest() {
        int longest_id = 0;
        for (int i = 0 ; i < All_seqs.size(); i++) {
            if (All_seqs[i].lens >= All_seqs[longest_id].lens) {
                longest_id = i;
            }
        }
        return longest_id;
    }


    // 直接并行计算全部k-mer-list的余弦距离
    //  余弦距离 ：
    //      respoint += (long) kmerA[i] * kmerB[i];
    //      lenA += (long) kmerA[i] * kmerA[i];
    //      lenB += (long) kmerB[i] * kmerB[i];
    //

    void Clu::get_cluster_two() {

        PLOGD << "======== this is Cluster::get_cluster_two ======";
        auto CPU_start = std::chrono::high_resolution_clock::now();

        tbb::parallel_for_each(All_seqs.begin() , All_seqs.end() , [&](Seq& seq) {
            long long tmp = 0;
            for (int i = 0 ;i < seq.k_mer_list.size() ; i++) {
                tmp += seq.k_mer_list[i] * seq.k_mer_list[i];
            }
            
            seq.A_A = tmp;
        });



        auto caculte_fcgr = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> fcgr_time = caculte_fcgr - CPU_start;
        PLOGD << "fcgr time: " << fcgr_time.count() << " ms";
        
        // 先按长度进行聚类
        int seq_count = All_seqs.size();
        std::map<int, std::vector<int>> Lens_clusters;
        
        while (seq_count) {
            int longest_id = Clu::selectLongest();
            Lens_clusters[longest_id] = Clu::cluster_for_lens(longest_id);
            PLOGD << "longest_id : " << longest_id;
            PLOGD << "Lens_clusters[id].size() : " << Lens_clusters[longest_id].size();
            seq_count -= Lens_clusters[longest_id].size();
        }

        PLOGD << "Lens_clusters.size() : "<<Lens_clusters.size() << std::endl;
        PLOGD << "======= START SECOND Cluster ======";
        // 对每个长度再进行聚类

        int tmp_count = 0;
        for (auto tmp : Lens_clusters) {

            PLOGD << "this lens_clusters count : " << tmp.second.size();

            // for(int i = 0; i < tmp.second.size(); i++){
            //     PLOGD << "seq id : " << tmp.second[i];
            // }


            FCGR_CU::one_cluster one_cluster(0.9,tmp.second);
            one_cluster.genClusters_CPU();

            for (auto tmp2 : one_cluster.cluster) {
                PLOGD << "cluster id : " << tmp2.first;
                PLOGD << "cluster count : " << tmp2.second.size();
                this->TOP_clusters[tmp2.first] = tmp2.second;
                tmp_count++;
            }
        }
        PLOGD << "======= FINISH Cluster =======";
        PLOGD << "cluster.size() : " << tmp_count;

        auto CPU_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> CPU_time = CPU_end - CPU_start;
        PLOGD << "CPU time: " << CPU_time.count() << " ms";
    }

    std::vector<int> get_seq_id(std::vector<std::string> &list){
        std::vector<int> res;
        for(auto &tmp : list){
            res.push_back(std::stoi(tmp));
        }
        return res;
    }
    void Clu::get_align_sort(){
        PLOGD << "==== out_seq_distance_matrix ====";
        std::vector<int> seq_id_list;
        for(auto &tmp : TOP_clusters){
            seq_id_list.push_back(tmp.first);
        }
        for(int i = 0 ; i < seq_id_list.size() ; i ++){
            for(int j = i + 1; j < seq_id_list.size() ; j ++){
                int col = seq_id_list[i];
                int row = seq_id_list[j];
                if(distance_matrix[col][row] == 0){
                    long long reP = FCGR_CU::get_Respoint(All_seqs[col].k_mer_list, All_seqs[row].k_mer_list);
                    long long tmpA = All_seqs[col].A_A;
                    long long tmpB = All_seqs[row].A_A;
                    double tmp_smi = reP / (sqrt(tmpB)*sqrt(tmpA));
                    distance_matrix[col][row] = 1 - tmp_smi;
                    distance_matrix[row][col] = 1 - tmp_smi;
                }
            }
        }
        std::string filename = "distance_matrix.txt";
        std::ofstream outfile(filename);
        if(!outfile.is_open()) {
            printf("open file failed");
            exit(0);
        }
        outfile << seq_id_list.size() << std::endl;
        for(int i = 0 ; i < seq_id_list.size() ; i ++){
            outfile <<std::left << std::setw(15) << seq_id_list[i];
            for(int j = 0; j < seq_id_list.size() ; j ++){
                outfile <<std::left << std::setw(10) << distance_matrix[seq_id_list[i]][seq_id_list[j]];
            }
            outfile << std::endl;
        }

        outfile.close();

        // ./decenttree -in ../example/example.dist -t NJ-R-V -no-banner -out njrv.newick
        std::string cmnd = "./include/decenttree";
        std::string out_res_file = "out.tree";
        cmnd.append(" ").append("-in ").append("distance_matrix.txt ").append("-t NJ-R-V ").append("-no-banner ").append("-out ").append(out_res_file);
        int res = system(cmnd.c_str());
        if (res != 0)
        {
            PLOGD << "run align softward false";
            exit(0);
        }else{
            PLOGD << "Tree file is ok: " << out_res_file;
        }
        std::ifstream infile(out_res_file);
        // ((2:0.066796586,(9:0.20109838,10:0.17624962):0.09366791):0.0047964305,(6:0.10352276,8:0.098947249):0.0064879432,7:0.11269307); // 生成一颗树
        std::string line = "";
        std::getline(infile, line);
        PLOGD << line;
        infile.close();

        std::map<int,std::vector<int>> process_list;
        int seq_flag = 0;
        std::vector<std::string> stack_str;  // 数据栈
        for (int i = 0; i < line.size(); i++)
        {
            if (line[i] == '(') // 左括号入栈
            {
                stack_str.push_back("(");
            }
            else if (line[i] == ')') // 右括号，弹出元素
            {
                std::vector<std::string> temp;
                while (stack_str.back() != "(")  // 一直弹出直到遇到左括号
                {
                    temp.push_back(stack_str.back());
                    stack_str.pop_back();
                }
                stack_str.pop_back();  // 弹出左括号

                seq_flag ++;
                // Cluster::get_seq_id(temp);
                process_list[seq_flag] = Cluster::get_seq_id(temp);
                // std::cout << "----------------" << std::endl;
            }
            else if (line[i] == ',') // 逗号分隔的元素
            {
                // 逗号直接跳过，处理后面的元素
                continue;
            }
            else // 处理元素（例如: '>_i10267:0.004044'）
            {
                int j = i;
                while (j < line.size() && line[j] != ',' && line[j] != ')')  // 直到遇到逗号或右括号
                {
                    j++;
                }
                // 将元素从i到j-1截取出来并压入栈
                std::string elem = line.substr(i, j - i);
                if(elem.find(":") == 0 || elem == ";"){
                    i = j - 1;  // 更新i，跳过已处理的部分
                    continue;
                }
                PLOGD << "elem : " << elem;
                stack_str.push_back(elem);  // 将元素压入栈

                i = j - 1;  // 更新i，跳过已处理的部分
            }
        }


        TOP_clusters_sort = process_list;
    }
}
