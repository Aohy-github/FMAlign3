#include "top.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, char* argv[]) {
    plog::init(plog::debug, "logfile.txt",100000, 1);
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

    if(classify == "false"){
        // 直接进行中心序列比对
        std::vector<Process_seq> process_seq_list;
        for(int i = 0; i < All_seqs.size(); i++){
            process_seq_list.push_back(Process_seq(i,All_seqs[i].content));
        }
        mode_pick(process_seq_list);

        for(int i = 0 ; i < process_seq_list.size() ; i ++){
            All_seqs[process_seq_list[i].seq_id].content = process_seq_list[i].content;
            // PLOGD << process_seq_list[i].content;
        }
    }else{
        // 执行中心+渐进式比对

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
                std::vector<Process_seq> process_seq_list;
                for(int i = 0; i < tmp.second.size(); i++){
                    process_seq_list.push_back(Process_seq(tmp.second[i],All_seqs[tmp.second[i]].content));
                }
                sort(process_seq_list.begin() ,  process_seq_list.end(),[](Process_seq &A , Process_seq &B){
                    return A.content.size() > B.content.size();
                });
                mode_pick(process_seq_list);
                for(int i = 0 ; i < process_seq_list.size() ; i ++){
                    All_seqs[process_seq_list[i].seq_id].content = process_seq_list[i].content;
                }
            }
        });

        PLOGD << "++ finish all align cluster ++";
        
        /**
         *  1.从每个群组中一条序列，计算剩余序列的距离矩阵 
         *  2.根据距离矩阵，计算指导树
         *  3.根据指导树进行序列比对（合并空格，对于新产生的空格，收集所有的空格情况，插回群组中）
         * 
         */
        cluster.get_align_sort();
        
        tbb::parallel_for_each(cluster.TOP_clusters.begin(), cluster.TOP_clusters.end(), [&](auto& tmp){
            if(tmp.second.size() > 1){
                std::vector<int> seq_id_list = tmp.second;
                std::string seq1 = All_seqs[seq_id_list[0]].content;
                Cluster::refresh_seq_content(seq_id_list,seq1);
            }
        });
    }
    
    out_sequence(outfilename);
    
   
    auto end_time = std::chrono::high_resolution_clock::now();



    std::chrono::duration<double, std::milli> top_end = end_time - top_time;
    PLOGD << "===== FINISH =====" << "time : " << top_end.count() << " ms";
    std::cout << "===== FINISH =====" << "time : " << top_end.count() << " ms\n";
    return 0;
}
