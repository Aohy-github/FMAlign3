#include "top.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;
void out_distance_matrix(Cluster::Clu& cluster){
    std::vector<std::vector<double>> new_distance_matrix;
    const int name_width = 10; // 序列名称列宽
    const int value_width = 12; // 数值列宽
    const int precision = 6; // 小数位数
    std::set<int> is_set;    // 已经处理的群组
    std::vector<int> matrix_id_list;
    int count = 0;
    
    for (auto& tmp : cluster.TOP_clusters) {
        if (tmp.second.size() > 1)
        {
            for(int i = 0; i < tmp.second.size(); i++){
                is_set.insert(tmp.second[i]);
                std::cout << tmp.second[i] << " ";
            }
            count += tmp.second.size() - 1;
        }
    }
    new_distance_matrix.resize(All_seqs.size() - count, std::vector<double>(All_seqs.size() - count, 0));

    //打印每一行
    int row = 0 , col = 0;
    std::cout << std::endl;
    std::vector<int> seq_id_list;
    for(int i = 0; i < distance_matrix.size(); i++){
        // 打印行标签，左对齐
        if(is_set.find(i)!= is_set.end()){
            seq_id_list.push_back(i);
            std::cout << std::left << std::setw(name_width) << All_seqs[i].name.substr(0,10) << " ";
            col = 0;
            
            for(int j = 0; j < distance_matrix[i].size(); j++){
                // 设置数值格式：固定小数位，右对齐
                if(is_set.find(j)!= is_set.end()){
                    if(distance_matrix[i][j] == 0 && i != j){
                        long long reP = FCGR_CU::get_Respoint(All_seqs[i].k_mer_list, All_seqs[j].k_mer_list);
                        long long tmpA = All_seqs[i].A_A;
                        long long tmpB = All_seqs[j].A_A;
                        double tmp_smi = reP / (sqrt(tmpB)*sqrt(tmpA));
                        distance_matrix[i][j] = 1 - tmp_smi;
                    }
                    new_distance_matrix[row][col++] = distance_matrix[i][j];

                    std::cout << std::right << std::fixed << std::setprecision(precision) 
                        << std::setw(value_width) << distance_matrix[i][j] << " ";
                }
                
            }
            std::cout << std::endl;
            row++;
        }
    }
    std::cout << std::endl;
    std::cout << std::endl;
    for(int i = 0; i < seq_id_list.size(); i++){
        // 打印行标签，左对齐
        
            std::cout << std::left << std::setw(name_width) << All_seqs[seq_id_list[i]].name.substr(0,10) << " ";
            col = 0;
            row++;
            for(int j = 0; j < new_distance_matrix[i].size(); j++){
                // 设置数值格式：固定小数位，右对齐
                    std::cout << std::right << std::fixed << std::setprecision(precision) 
                        << std::setw(value_width) << new_distance_matrix[i][j] << " ";
                
            }
            std::cout << std::endl;
    }

    //打印表头
    std::cout << std::endl;
    std::cout << std::setw(name_width) << " " << " "; // 左上角空白
    for(int i = 0; i < distance_matrix.size(); i++){
        if(is_set.find(i) != is_set.end())
        std::cout << std::setw(value_width) << All_seqs[i].name.substr(0,10) << " ";
    }
        
    std::cout << std::endl;
}

int main(int argc, char* argv[]) {
    plog::init(plog::debug, "logfile.txt" , 100000, 1);
    auto top_time = std::chrono::high_resolution_clock::now();
    

    int kmer_size;
    std::string classify;
    std::string filename;
    std::string outfilename;
    // data/truncated_sequences_1500.fasta
    // data/Ecoli.fasta
    //steps-512_0_NOT_U.fasta
    //23s_rRNA.fasta
    // chr1.fasta
    
    
    // //std::string out_filename = "../data/FCGR.txt";
    
    po::options_description desc("选项:");
    desc.add_options()
        ("help,h", "显示帮助信息")
        ("input,i", po::value<std::string>(&filename)->required(), "输入 FASTA 文件路径")
        ("kmer,k", po::value<int>(&kmer_size)->default_value(15), "k-mer 大小")
        ("classify,c", po::value<std::string>(&classify)->default_value("false"), "是否进行分类")
        ("output,o", po::value<std::string>(&outfilename)->default_value("output.fasta"), "输出文件路径");

    po::variables_map vm;
    try {
        // 解析命令行参数
        po::store(po::parse_command_line(argc, argv, desc), vm);

        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 0;
        }

        // 检查必需参数
        po::notify(vm);
    }
    catch(po::error& e) {
        std::cerr << "错误: " << e.what() << "\n";
        std::cerr << desc << "\n";
        return 1;
    }
    
    

    PLOGD << "k-mer : " << kmer_size;
    PLOGD << "classify : " << classify;
    PLOGD << "input : " << filename;
    PLOGD << "output : " << outfilename;
    PLOGD << "";
    PLOGD << "";

    // 获取序列
    All_seqs = get_sequence(filename);
    PLOGD << "get seqs : "<<All_seqs.size();

    if(classify == "false"){
        // 直接进行中心序列比对
        std::vector<int> seq_id_list;
        for(int i = 0; i < All_seqs.size(); i++){
            seq_id_list.push_back(i);
        }
        mode_pick(seq_id_list);
    }else{
        // 执行中心+渐进式比对
        PLOGD << "===== START FCGR =====";

        auto CPU_start = std::chrono::high_resolution_clock::now();
        oneapi::tbb::parallel_for_each(All_seqs.begin(), All_seqs.end(), [&](Seq& one){
            bool bo = FCGR_CU::get_one_FCGR_CPU(one);
            if(!bo){
                std::cout << one.name << "false !" << std::endl;
            }
            //PLOGD<<"kmer.size()  : " << one.k_mer_list.size();
        });

        auto CPU_end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> CPU_time = CPU_end - CPU_start;

        
        // // exit(0);
        PLOGD << "===== FINISH FCGR =====" << "CPU time : " << CPU_time.count() << " ms";
        PLOGD << "";
        PLOGD << "";
        PLOGD <<std::setw(10)<< std::setfill(' ') <<"   =====   START CLUSTER   =====   ";
        PLOGD << "";
        PLOGD << "";
        auto cluster_start = std::chrono::high_resolution_clock::now();
        
        Cluster::Clu cluster;
        cluster.get_cluster_two();

        auto cluster_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> cluster_time = cluster_end - cluster_start;
        

        // int count = 0;
        // for (auto tmp : cluster.TOP_clusters) {
        //     if (tmp.second.size() > 1) {
        //         std::cout <<" size > 2 : "<< tmp.second.size() << std::endl;
        //         count += tmp.second.size();
        //     }
        // }
        // std::cout << "count : "<< count << std::endl;
        PLOGD << "";
        PLOGD << "";
        PLOGD << "cluster size : " << cluster.TOP_clusters.size();
        PLOGD << "";
        PLOGD << "";
        PLOGD << "===== FINISH CLUSTER =====" << "time : " << cluster_time.count() << " ms";
        PLOGD << "";
        PLOGD << "";

        tbb::parallel_for_each(cluster.TOP_clusters.begin(), cluster.TOP_clusters.end(), [&](auto& tmp){
            if(tmp.second.size() > 1){
                std::vector<int> seq_id_list = tmp.second;
                sort(tmp.second.begin() ,  tmp.second.end());
                mode_pick(tmp.second);
            }
        });

        PLOGD << "++ finish all align cluster ++";
        
        // 整理所有的群组
        cluster.get_align_sort();

        for(auto tmp : cluster.TOP_clusters_sort){
            PLOGD << "list id : " << tmp.first;
            PLOGD << "cluster count : " << tmp.second.size();
            for(int i = 0 ; i < tmp.second.size(); i++){
                PLOGD << "seq id : " << tmp.second[i];
            }
        }
        // out_distance_matrix(cluster);
        // 生成新的距离矩阵
        /**
         *  1.从每个群组中一条序列，计算剩余序列的距离矩阵 
         *  2.根据距离矩阵，计算指导树
         *  3.根据指导树进行序列比对（合并空格，对于新产生的空格，收集所有的空格情况，插回群组中）
         * 
         */
        
    }
    
    out_sequence(outfilename);
    
   
    auto end_time = std::chrono::high_resolution_clock::now();



    std::chrono::duration<double, std::milli> top_end = end_time - top_time;
    PLOGD << "===== FINISH =====" << "time : " << top_end.count() << " ms";
    std::cout << "===== FINISH =====" << "time : " << top_end.count() << " ms\n";
    return 0;
}
