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
        
        
        // 先按长度进行聚类
        int seq_count = All_seqs.size();
        std::map<int, std::vector<int>> Lens_clusters;
        
        while (seq_count) {
            int longest_id = Clu::selectLongest();
            Lens_clusters[longest_id] = Clu::cluster_for_lens(longest_id);
            // PLOGD << "longest_id : " << longest_id;
            // PLOGD << "Lens_clusters[id].size() : " << Lens_clusters[longest_id].size();
            seq_count -= Lens_clusters[longest_id].size();
        }

        // 对每个长度再进行聚类

        int tmp_count = 0;
        for (auto tmp : Lens_clusters) {
            // for(int i = 0; i < tmp.second.size(); i++){
            //     PLOGD << "seq id : " << tmp.second[i];
            // }


            FCGR_CU::one_cluster one_cluster(0.9,tmp.second);
            one_cluster.genClusters_CPU();

            for (auto tmp2 : one_cluster.cluster) {
                // PLOGD << "cluster id : " << tmp2.first;
                // PLOGD << "cluster count : " << tmp2.second.size();
                this->TOP_clusters[tmp2.first] = tmp2.second;
                tmp_count++;
            }
        }


        auto CPU_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> CPU_time = CPU_end - CPU_start;

    }


    
    std::vector<int> get_seq_id(std::vector<std::string> &list){
        std::vector<int> res;
        for(auto &tmp : list){
            res.push_back(std::stoi(tmp));
        }
        return res;
    }



    std::string get_common_seq(std::vector<std::string> content_list){
        
        size_t seq_length = content_list[0].length();
        std::string final_seq(seq_length , 'N');
        std::vector<int> N_option;
        for(int i = 0; i < content_list.size() ; i++){
            if(content_list[i].length() != seq_length){
 
                exit(0);
            }
        }
        for(size_t i = 0; i < seq_length ; ++i){
            std::unordered_map<char,int> char_map;
            for(auto &tmp : content_list){
                char_map[tmp[i]]++;
            }
            char max_char = ' ';
            int max_count = 0;
            for(auto &tmp : char_map){
                if(tmp.second > max_count){
                    max_count = tmp.second;
                    max_char = tmp.first;
                }
            }
            if(max_char != '-' && max_count >= content_list.size()* 0.6){
                final_seq[i] = max_char;
            }else {
                final_seq[i] = (content_list[0][0] >= 'A') ? 'N' : 'n';
                N_option.push_back(i);
            }
        }
        // for(int i = 0; i < N_option.size() - 1; i++){
        //     int index = N_option[i];
        //     int next = N_option[i + 1];
        //     if(next - index < 21){
        //         for(int j = index; j < next; j++){
        //             final_seq[i] = (content_list[0][0] >= 'A') ? 'N' : 'n';
        //         }
        //     }
        // }
        
        return final_seq;
    }

    void profile_more_seq(std::vector<int> seq_id_list){
        std::vector<Process_seq> process_seq_list;
        for(int i = 0 ;i < seq_id_list.size() ; i++){
            process_seq_list.push_back(Process_seq(seq_id_list[i],All_seqs[seq_id_list[i]].content));
        }
        std::sort(process_seq_list.begin() , process_seq_list.end() , [](Process_seq &A , Process_seq &B){
            return A.content.size() > B.content.size();
        });
        
        std::vector<std::vector<Align_Seq>> not_split_seq_list;
        std::vector<Align_Seq> tmp;
        for(int i = 0; i < process_seq_list.size(); i++){
            tmp.push_back(Align_Seq(process_seq_list[i].content));
            process_seq_list[i].content.clear();
        }
        not_split_seq_list.push_back(tmp);
        star_align(not_split_seq_list);
        merge_split(not_split_seq_list,process_seq_list);

        for(int i = 0 ; i < process_seq_list.size() ; i++){
            All_seqs[process_seq_list[i].seq_id].content = process_seq_list[i].content;
        }
        
    }
    void refresh_seq_content(std::vector<int>& seq_id_list , std::string seq1){
        tbb::parallel_for_each(seq_id_list.begin(), seq_id_list.end(), [&](int id){
            std::string tmp = All_seqs[id].content;
            std::string final_res(seq1.size() , '-');
            int le = 0, ri = 0;
            while(le < seq1.size() && ri < tmp.size()){
                if(seq1[le] == '-' && tmp[ri] != '-'){ // le != ri
                    final_res[le] = '-';
                    le++;
                }else{
                    final_res[le] = tmp[ri];
                    ri++;
                    le++;
                }
            }
            All_seqs[id].content = final_res;
        });
    }
    void Clu::get_align_sort(){
        std::vector<int> process_id_list;
        for(auto &tmp : TOP_clusters){
            process_id_list.push_back(tmp.first);
            //PLOGD << "tmp.first :" << tmp.first << " seq name :" << All_seqs[tmp.first].name;
        }
        for(int i = 0 ; i < process_id_list.size() ; i ++){
            for(int j = i + 1; j < process_id_list.size() ; j ++){
                int col = process_id_list[i];
                int row = process_id_list[j];
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
        outfile << process_id_list.size() << std::endl;
        for(int i = 0 ; i < process_id_list.size() ; i ++){
            outfile <<std::left << std::setw(15) << process_id_list[i];
            for(int j = 0; j < process_id_list.size() ; j ++){
                outfile <<std::left << std::setw(10) << distance_matrix[process_id_list[i]][process_id_list[j]];
            }
            outfile << std::endl;
        }

        outfile.close();

        // ./decenttree -in ../example/example.dist -t NJ-R-V -no-banner -out njrv.newick
        std::string cmnd = "./include/decenttree ";
        std::string out_res_file = "out.tree";
        cmnd.append("-in ").append("distance_matrix.txt ")
            .append("-t BIONJ-R ")
            .append("-no-banner ")
            .append("-out ").append(out_res_file).append("> /dev/null");

        int res = system(cmnd.c_str());
        if (res != 0){
            PLOGD << "run align softward false";
            exit(0);
        }else{
            PLOGD << "Tree file is ok: " << out_res_file;
        }
        std::ifstream infile(out_res_file);
        
        std::string line = "";
        std::getline(infile, line);
        infile.close();

        std::vector<std::vector<int>> mid_seq_id;
        
        std::vector<std::string> stack_str;  

        for (int i = 0; i < line.size(); i++)
        {
            if (line[i] == '(') {
                stack_str.push_back("(");
            }
            else if (line[i] == ')'){
                std::vector<std::string> temp;
                while (stack_str.back() != "(") // 如果弹出的元素数量是大于2个，那就说明是一个节点，进行比对，否则就从mid_seq_content中弹出一个元素进行比对
                {
                    temp.push_back(stack_str.back());
                    stack_str.pop_back();
                }
                stack_str.pop_back();  // 弹出左括号
                
                if(temp.size() >= 2){
                    
                    std::vector<int> seq_id_list = Cluster::get_seq_id(temp);
                    Cluster::profile_more_seq(seq_id_list);
                    mid_seq_id.push_back(seq_id_list);

                }else if(temp.size() == 1){
                   
                    std::vector<std::string> content_list;
                    std::vector<int> seq_id_list = mid_seq_id.back();
                    mid_seq_id.pop_back();
                    
                    for(int j = 0 ; j < seq_id_list.size() ; j++){
                        content_list.push_back(All_seqs[seq_id_list[j]].content);
                        
                    }
                    std::string seq1 = Cluster::get_common_seq(content_list);
                    std::string seq2 = All_seqs[std::stoi(temp[0])].content;
                    
                    align_wfa_profile(seq2 , seq1);
                    
                    Cluster::refresh_seq_content(seq_id_list,seq1);
                    
                    All_seqs[std::stoi(temp[0])].content = seq2;
                    seq_id_list.push_back(std::stoi(temp[0]));
                    mid_seq_id.push_back(seq_id_list);

   
                }else{
                    
                    std::vector<int> seq_id_list1 = mid_seq_id.back();
                    mid_seq_id.pop_back();
                    std::vector<int> seq_id_list2 = mid_seq_id.back();
                    mid_seq_id.pop_back();

                    std::vector<std::string> content_list1;
                    for(int j = 0; j < seq_id_list1.size() ; j++){
                        content_list1.push_back(All_seqs[seq_id_list1[j]].content);
                    }
                    std::string seq1 = Cluster::get_common_seq(content_list1);
                    int len1 = seq1.size();
                    
                    std::vector<std::string> content_list2;
                    
                    for(int j = 0; j < seq_id_list2.size() ; j++){
                        content_list2.push_back(All_seqs[seq_id_list2[j]].content);
                        
                    }
                    std::string seq2 = Cluster::get_common_seq(content_list2);
                    int len2 = seq2.size();
                    
                    align_wfa_profile(seq2 , seq1);

                    if(len1 != seq1.size()){
                        Cluster::refresh_seq_content(seq_id_list1,seq1);
                    }
                    if(len2 != seq2.size()){
                        Cluster::refresh_seq_content(seq_id_list2,seq2);
                    }
                    
                    seq_id_list1.insert(seq_id_list1.end(),seq_id_list2.begin(),seq_id_list2.end());
                    mid_seq_id.push_back(seq_id_list1);
 
                }
            }
            else if (line[i] == ',') {// 逗号分隔的元素
                continue;
            }
            else {
                int j = i;
                while (j < line.size() && line[j] != ',' && line[j] != ')' &&line[j] != '(')  // 直到遇到逗号或右括号
                {
                    j++;
                }
                // 将元素从i到j-1截取出来并压入栈
                std::string elem = line.substr(i, j - i);
                
                if(elem.find(":") == 0 || elem == ";"){
                    i = j - 1;  // 更新i，跳过已处理的部分
                    continue;
                }
                stack_str.push_back(elem);  // 将元素压入栈
                
                i = j - 1;  // 更新i，跳过已处理的部分
            }
        }
        std::vector<int> seq_id_list1 = mid_seq_id.back();
        mid_seq_id.pop_back();

        std::vector<std::string> content_list1;
        for(int j = 0; j < seq_id_list1.size() ; j++){
            content_list1.push_back(All_seqs[seq_id_list1[j]].content);
            
        }
        std::string seq1 = Cluster::get_common_seq(content_list1);
        int len1 = seq1.size();

        std::vector<int> seq_id_list2 = mid_seq_id.back();
        mid_seq_id.pop_back();
        std::vector<std::string> content_list2;
        for(int j = 0; j < seq_id_list2.size() ; j++){
            content_list2.push_back(All_seqs[seq_id_list2[j]].content);
        }
        std::string seq2 = Cluster::get_common_seq(content_list2);
        int len2 = seq2.size();
        align_wfa_profile(seq2 , seq1);
        
        if(seq1.size() != len1){
            Cluster::refresh_seq_content(seq_id_list1,seq1);
        }
        if(seq2.size() != len2){
            Cluster::refresh_seq_content(seq_id_list2,seq2);
        }
    }
   
}


