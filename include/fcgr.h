//
// Created by aohy on 24-12-24.
//

#ifndef FCGR_CUH
#define FCGR_CUH

#include "top.h"


// 并行计算FCGR

namespace FCGR_CU {
    // bool get_one_FCGR_GPU(Seq& seq);
    bool get_one_FCGR_CPU(Seq& seq);
    long long get_Respoint(std::vector<int>& A_kmer_list, std::vector<int>& B_kmer_list);
    // __global__ void computeProducts(const int* kmerA, const int* kmerB, double* dotProduct, double* sumA, double* sumB, int n);
    // __global__ void cuda_processSeqFCGR(int n, const char* str, int K, int* k_mer_Matrix);
    class one_cluster {
    public:
        double sim;
        std::vector<int> total_SeqList;
        std::unordered_map<int,bool> hashMap;
        std::map<int, std::vector<int>> cluster;
        one_cluster(double sim , std::vector<int>& id_list) : sim(sim) , total_SeqList(id_list) {};
        // cuda计算距离
        //double getDistanceCUDA(const std::vector<int>& kmerA, const std::vector<int>& kmerB);
        // 对当前seq进行聚类GPU
        //void genClusters_GPU();
        // 对当前seq进行聚类CPU
        void genClusters_CPU();
        // 获取最长序列id
        int getLoneID();
        // 根据long_id进行聚类
        //std::vector<int> get_smi_seqs_GPU(int long_id);
        std::vector<int> get_smi_seqs_CPU(int long_id);
    };
}


// #define CUDA_CHECK(call)                                                     \
//     do {                                                                     \
//         cudaError_t error = call;                                            \
//         if (error != cudaSuccess) {                                          \
//             std::cerr << "CUDA Error: " << cudaGetErrorString(error)         \
//                       << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
//             exit(EXIT_FAILURE);                                              \
//         }                                                                    \
// } while (0)



#endif //FCGR_CUH
