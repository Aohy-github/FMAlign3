#include "top.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, char* argv[]) {

    auto top_time = std::chrono::high_resolution_clock::now();
    

    int kmer_size;
    bool classify;
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
        ("classify,c", po::bool_switch(&classify)->default_value(false), "是否进行分类")
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

    plog::init(plog::debug, "logfile.txt"); // 设置日志级别和日志文件




    All_seqs = get_sequence(filename);
    // PLOGD << "get seqs : "<<All_seqs.size();



    // PLOGD << "===== START FCGR =====";

    // auto CPU_start = std::chrono::high_resolution_clock::now();
    // oneapi::tbb::parallel_for_each(All_seqs.begin(), All_seqs.end(), [&](Seq& one){
    //     bool bo = FCGR_CU::get_one_FCGR_CPU(one);
    //     if(!bo){
    //         std::cout << one.name << "false !" << std::endl;
    //     }
    //     //PLOGD<<"kmer.size()  : " << one.k_mer_list.size();
    // });

    // auto CPU_end = std::chrono::high_resolution_clock::now();

    // std::chrono::duration<double, std::milli> CPU_time = CPU_end - CPU_start;

    
    // // // exit(0);
    // PLOGD << "===== FINISH FCGR =====" << "CPU time : " << CPU_time.count() << " ms";
    // std::cout << "===== FINISH FCGR =====" << "CPU time : " << CPU_time.count() << " ms\n";
    // PLOGD << "";
    // PLOGD << "";
    // PLOGD <<std::setw(10)<< std::setfill(' ') <<"   =====   START CLUSTER   =====   ";
    // PLOGD << "";
    // PLOGD << "";
    // auto cluster_start = std::chrono::high_resolution_clock::now();
    
    
    
    // Cluster::Clu cluster;
    // cluster.get_cluster_two();

    // auto cluster_end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double, std::milli> cluster_time = cluster_end - cluster_start;
    

    // // int count = 0;
    // // for (auto tmp : cluster.TOP_clusters) {
    // //     if (tmp.second.size() > 1) {
    // //         std::cout <<" size > 2 : "<< tmp.second.size() << std::endl;
    // //         count += tmp.second.size();
    // //     }
    // // }
    // // std::cout << "count : "<< count << std::endl;
    // PLOGD << "";
    // PLOGD << "";
    // PLOGD << "cluster size : " << cluster.TOP_clusters.size();
    // PLOGD << "";
    // PLOGD << "";
    // PLOGD << "===== FINISH CLUSTER =====" << "time : " << cluster_time.count() << " ms";
    // std::cout << "===== FINISH CLUSTER =====" << "time : " << cluster_time.count() << " ms\n";
    // PLOGD << "";
    // PLOGD << "";

    // PLOGD << "===== START ALIGN =====";
    // 进行比对 1.对每个群组进行比对，
    // for(auto& tmp : cluster.TOP_clusters) {
    //     PLOGD << "tmp.first : " << tmp.first;
    //     PLOGD << "tmp.second.size() : " << tmp.second.size();
    //     sort(tmp.second.begin() ,  tmp.second.end()); // 整理编号顺序
        
    //     if(tmp.second.size() > 1){
    //         for(int i = 0; i < tmp.second.size(); i++){
    //             PLOGD << tmp.second[i] << " seq_name : " << All_seqs[tmp.second[i]].name << "lens : " << All_seqs[tmp.second[i]].content.size();
    //         }
    //         center_align(cluster.TOP_clusters[tmp.first]);
    //     }
            
    // }
    // bool bo = false;
    
    // if(bo){
    //     tbb::parallel_for_each(cluster.TOP_clusters.begin(), cluster.TOP_clusters.end(), [&](auto& tmp){
    //         std::vector<int> seq_id_list = tmp.second;
    //         sort(tmp.second.begin() ,  tmp.second.end());
    //         if(tmp.second.size() > 1)
    //             center_align(tmp.second);
    //     });
    // }else{
    //     for(auto& tmp :cluster.TOP_clusters){
    //     PLOGD << "tmp.first : " << tmp.first;
    //     PLOGD << "tmp.second.size() : " << tmp.second.size();
    //     std::vector<int> seq_id_list = tmp.second;
    //     if(tmp.second.size() > 1){
    //         for(int i = 0; i < tmp.second.size(); i++){
    //             PLOGD << tmp.second[i] << " seq_name : " << All_seqs[tmp.second[i]].name << "lens : " << All_seqs[tmp.second[i]].content.size();
    //         }
    //         center_align(tmp.second);
    //     }
    // }
    // }
    // // 处理每个群组
    // for(int i = 0; i < distance_matrix.size(); i++)
    //     std::cout <<std::setw(15)<< All_seqs[i].name.substr(0,5) << " ";
    // std::cout << std::endl;
    // for(int i = 0; i < distance_matrix.size(); i++){
    //     std::cout<< All_seqs[i].name.substr(0,5) << " ";
    //     for(int j = 0; j < distance_matrix.size(); j++){
    //         std::cout <<std::setw(13)<< distance_matrix[i][j];
    //     }
    //     std::cout << std::endl;
    // }
    
    std::vector<int> seq_id_list;
    for(int i = 0; i < All_seqs.size(); i++){
        seq_id_list.push_back(i);
    }
    process_more_seq_and_lens_less_10000(seq_id_list);
    out_sequence(outfilename);
    
   
    auto end_time = std::chrono::high_resolution_clock::now();



    std::chrono::duration<double, std::milli> top_end = end_time - top_time;
    PLOGD << "===== FINISH =====" << "time : " << top_end.count() << " ms";
    std::cout << "===== FINISH =====" << "time : " << top_end.count() << " ms\n";
    return 0;
}
