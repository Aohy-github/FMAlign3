//
// Created by aohy on 24-12-24.
//

#include "fcgr.h"

namespace FCGR_CU {
    bool get_one_FCGR_CPU(Seq& seq) {
        int k = 6;
        int mask = (1 << (2 * k)) - 1;
        std::array<int, 256> char_to_code = {};
        char_to_code['A'] = 0;
        char_to_code['T'] = 1;
        char_to_code['C'] = 2;
        char_to_code['G'] = 3;

        char_to_code['a'] = 0;
        char_to_code['t'] = 1;
        char_to_code['c'] = 2;
        char_to_code['g'] = 3;
        std::vector<int> k_mer_CPU(seq.k_mer_list.size() , 0);
        std::string content = seq.content;
        int lens = content.size();
        unsigned int hash_val = 0;

        for (int i = 0; i < k; i++) { //读取前 k 位碱基
            hash_val = (hash_val << 2) | char_to_code[content[i]];
        }

        hash_val &= mask;
        k_mer_CPU[hash_val & mask]++;
        
        //std::cout << std::endl;
        for (int j = 6; j < lens; j++) {
            hash_val = (hash_val << 2) | char_to_code[content[j]];;
            hash_val &= mask;
            k_mer_CPU[hash_val]++;
        }
        seq.k_mer_list = std::vector<int>(k_mer_CPU.begin(), k_mer_CPU.end());
        return true;
    }
 
    void one_cluster::genClusters_CPU() {
        for (int i = 0 ;i < this->total_SeqList.size() ; i++) {
            // PLOGD << "this->total_SeqList : " << this->total_SeqList[i];
            this->hashMap[this->total_SeqList[i]] = false;
        }
        int seq_count = this->total_SeqList.size();

        while (seq_count >= 1) {
            int long_id = getLoneID();
            // PLOGD << "long_id : " << long_id;
            if (long_id >= 0) {
                this->cluster[long_id] = get_smi_seqs_CPU(long_id);
                seq_count -= this->cluster[long_id].size();
            }else {
                PLOGD << "==== 数组越界！====";
                exit(0);
            }
        }

    }

    std::vector<int> one_cluster::get_smi_seqs_CPU(int long_id) {
        std::vector<int> res;
        res.push_back(long_id);
        this->hashMap[long_id] = true;
        long long tmpA = All_seqs[long_id].A_A;
        for (int i = 0; i < total_SeqList.size(); i++) { // 从 0 开始会包含 long_id 自身

            if (!this->hashMap[this->total_SeqList[i]]) {
                long long tmpB = All_seqs[this->total_SeqList[i]].A_A;
                long long reP = get_Respoint(All_seqs[long_id].k_mer_list , All_seqs[this->total_SeqList[i]].k_mer_list);
                double tmp_smi = reP / (sqrt(tmpB)*sqrt(tmpA));
                //std::cout << "tmp_smi : " << tmp_smi << std::endl;
                
                distance_matrix[long_id][this->total_SeqList[i]] = 1 - tmp_smi;
                distance_matrix[this->total_SeqList[i]][long_id] = 1 - tmp_smi;

                if (tmp_smi >= this->sim) {
                    res.push_back(this->total_SeqList[i]);
                    this->hashMap[this->total_SeqList[i]] = true;
                }
            }
        }
        return res;
    }

    long long get_Respoint(std::vector<int>& A_kmer_list, std::vector<int>& B_kmer_list) {
        long long res = 0 ;
        int lens = A_kmer_list.size();
        tbb::enumerable_thread_specific<long long> local_resPoints(0);
        // 4096 / 4 = 1024
        tbb::parallel_for(tbb::blocked_range<size_t>(0 , lens, 512),[&](const tbb::blocked_range<size_t>& r) {
            long long& local_resP = local_resPoints.local();
            for (int i = r.begin() ; i < r.end() ; i++) {
                local_resP += (A_kmer_list[i] * B_kmer_list[i]);
            }
        });
        for (auto tmp : local_resPoints) {
            res += tmp;
        }
        // PLOGD << "res : " << res << std::endl;
        return res;
    }


    int one_cluster::getLoneID() {
        int idx = -1 , maxLen = 0;
        for (int i = 0 ; i < this->total_SeqList.size() ; i++) {
            int id = this->total_SeqList[i];
            if (!this->hashMap[id] && All_seqs[id].content.size() > maxLen) {
                maxLen = All_seqs[id].content.size();
                idx = id;
            }
        }
        //std::cout << "long id: "<< idx << std::endl;
        return idx;
    }
}



/**
    auto GPU_start = std::chrono::high_resolution_clock::now();
    oneapi::tbb::parallel_for_each(All_seqs.begin(), All_seqs.end(), [&](Seq& one){
        bool bo = FCGR_CU::get_one_FCGR_GPU(one);
        if(!bo){
            std::cout << one.name << "false !" << std::endl;
        }
    });
    auto GPU_end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> gpu_time = GPU_end - GPU_start;

    std::cout << "GPU time: " << gpu_time.count() << " ms" << std::endl;
**/


   // bool get_one_FCGR_GPU(Seq& seq) {
    //     //创建GPU内存
    //     char* d_sequence;
    //     int* d_kmer_Matrix;
    //     int slens = seq.content.size();
    //     int num_kmers = 1 << (2*K); // 4^K

    //     cudaMalloc(&d_sequence, slens * sizeof(char));
    //     cudaMemcpy(d_sequence , seq.content.c_str() , slens * sizeof(char), cudaMemcpyHostToDevice);

    //     cudaMalloc(&d_kmer_Matrix, num_kmers * sizeof(int));
    //     cudaMemset(d_kmer_Matrix, 0, num_kmers * sizeof(int));

    //     int threads_per_block = 256; // 一共256个线程
    //     int blocks_per_grid = (slens + threads_per_block - 1) / threads_per_block; // slens(1596) + 256 - 1 / 256 大约 7 个块儿 ， 没块儿256个线程
    //     size_t shared_mem_size = num_kmers * sizeof(int); // 动态共享内存大小

    //     FCGR_CU::cuda_processSeqFCGR<<<blocks_per_grid, threads_per_block , shared_mem_size>>>(slens , d_sequence , K , d_kmer_Matrix);

    //     cudaDeviceSynchronize();
    //     cudaError_t err = cudaGetLastError();
    //     if(err != cudaSuccess) {
    //         std::cout << "Error: " << cudaGetErrorString(err) << std::endl;
    //         cudaFree(d_sequence);
    //         cudaFree(d_kmer_Matrix);
    //         exit(1);
    //     }

    //     cudaMemcpy(seq.k_mer_list.data() , d_kmer_Matrix, 4096 * sizeof(int), cudaMemcpyDeviceToHost);

    //     cudaFree(d_sequence);
    //     cudaFree(d_kmer_Matrix);
    //     return true;
    // }

    // double one_cluster::getDistanceCUDA(const std::vector<int>& kmerA, const std::vector<int>& kmerB) {
    //     if (kmerA.size() != kmerB.size()) {
    //         std::cerr << "Error: kmerA and kmerB must have the same length." << std::endl;
    //         return -1.0;
    //     }

    //     int n = kmerA.size();
    //     size_t size = n * sizeof(int);


    //     // 分配设备内存
    //     int* d_kmerA = nullptr;
    //     int* d_kmerB = nullptr;
    //     double* d_dotProduct = nullptr;
    //     double* d_sumA = nullptr;
    //     double* d_sumB = nullptr;

    //     // 分配设备内存
    //     CUDA_CHECK(cudaMalloc((void**)&d_kmerA, size));
    //     CUDA_CHECK(cudaMalloc((void**)&d_kmerB, size));
    //     CUDA_CHECK(cudaMalloc((void**)&d_dotProduct, sizeof(double)));
    //     CUDA_CHECK(cudaMalloc((void**)&d_sumA, sizeof(double)));
    //     CUDA_CHECK(cudaMalloc((void**)&d_sumB, sizeof(double)));

    //     // 初始化归约结果为 0
    //     CUDA_CHECK(cudaMemset(d_dotProduct, 0, sizeof(double)));
    //     CUDA_CHECK(cudaMemset(d_sumA, 0, sizeof(double)));
    //     CUDA_CHECK(cudaMemset(d_sumB, 0, sizeof(double)));

    //     // 复制数据到设备
    //     CUDA_CHECK(cudaMemcpy(d_kmerA, kmerA.data(), size, cudaMemcpyHostToDevice));
    //     CUDA_CHECK(cudaMemcpy(d_kmerB, kmerB.data(), size, cudaMemcpyHostToDevice));

    //     // 声明线程数量
    //     int threadsPerBlock = 256;
    //     //声明Grid中块儿的数量
    //     int blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;
    //     // int blocksPerGrid = 2;
    //     // 计算共享内存大小：3 * threadsPerBlock * sizeof(double)
    //     size_t sharedMemSize = 3 * threadsPerBlock * sizeof(double);

    //     // 启动核函数
    //     FCGR_CU::computeProducts<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(d_kmerA, d_kmerB, d_dotProduct, d_sumA, d_sumB, n);
    //     CUDA_CHECK(cudaDeviceSynchronize());

    //     // 检查核函数执行是否有错误
    //     CUDA_CHECK(cudaGetLastError());

    //     // 将结果从设备复制回主机
    //     double h_dotProduct = 0.0;
    //     double h_sumA = 0.0;
    //     double h_sumB = 0.0;

    //     CUDA_CHECK(cudaMemcpy(&h_dotProduct, d_dotProduct, sizeof(double), cudaMemcpyDeviceToHost));
    //     CUDA_CHECK(cudaMemcpy(&h_sumA, d_sumA, sizeof(double), cudaMemcpyDeviceToHost));
    //     CUDA_CHECK(cudaMemcpy(&h_sumB, d_sumB, sizeof(double), cudaMemcpyDeviceToHost));

    //     // 释放设备内存
    //     CUDA_CHECK(cudaFree(d_kmerA));
    //     CUDA_CHECK(cudaFree(d_kmerB));
    //     CUDA_CHECK(cudaFree(d_dotProduct));
    //     CUDA_CHECK(cudaFree(d_sumA));
    //     CUDA_CHECK(cudaFree(d_sumB));

    //     // 计算余弦相似度
    //     double sim = h_dotProduct / std::sqrt(h_sumA * h_sumB);
    //     return sim;
    // }
    // __global__ void computeProducts(const int* kmerA, const int* kmerB, double* dotProduct, double* sumA, double* sumB, int n) {
    //     extern __shared__ double sharedData[]; // 动态共享内存 每个 Thread block 中的所有 Thread 共用同一块 Shared memory

    //     unsigned int tid = threadIdx.x;
    //     unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;

    //     double tempDot = 0.0;
    //     double tempA = 0.0;
    //     double tempB = 0.0;

    //     if (index < n) {
    //         double valA = static_cast<double>(kmerA[index]);
    //         double valB = static_cast<double>(kmerB[index]);

    //         tempDot = valA * valB;  // A[i]*B[i]
    //         tempA   = valA * valA;  // A[i]*A[i]
    //         tempB   = valB * valB;  // B[i]*B[i]
    //     }else {
    //         tempDot = 0.0;
    //         tempA = 0.0;
    //         tempB = 0.0;
    //     }

    //     // 每个线程将结果存储在共享内存中
    //     sharedData[tid] = tempDot;
    //     sharedData[tid + blockDim.x] = tempA;
    //     sharedData[tid + 2 * blockDim.x] = tempB;
    //     __syncthreads();

    //     // 归约操作（块内求和）
    //     for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    //         if (tid < s) {
    //             sharedData[tid]                  += sharedData[tid + s];
    //             sharedData[tid + blockDim.x]     += sharedData[tid + blockDim.x + s];
    //             sharedData[tid + 2 * blockDim.x] += sharedData[tid + 2 * blockDim.x + s];
    //         }
    //         __syncthreads();
    //     }

    //     // 将每个块的部分和写回全局内存
    //     if (tid == 0) {
    //         atomicAdd(dotProduct, sharedData[0]);
    //         atomicAdd(sumA,       sharedData[blockDim.x]);
    //         atomicAdd(sumB,       sharedData[2 * blockDim.x]);
    //     }
    // }

    // __global__ void cuda_processSeqFCGR(int n ,const char* str, int K , int* k_mer_Matrix) {
    //     int num_kmers = 1 << (2*K); // 4^K = 2^(2K)
    //     n = n - K;
    //     int tid = blockIdx.x * blockDim.x + threadIdx.x; // 全局线程索引
    //     int stride = blockDim.x * gridDim.x;            // 跨步大小

    //     extern __shared__ int shared_counts[]; // 动态共享内存
    //     for (int i = threadIdx.x; i < num_kmers; i += blockDim.x) {
    //         shared_counts[i] = 0; // 初始化共享内存
    //     }
    //     __syncthreads();

    //         for (int i = tid; i <= n; i += stride) {
    //             int hash_val = 0;
    //             for (int j = 0; j < K; j++) {
    //                 char ch = str[i + j];
    //                 switch (ch) {
    //                     case 'A': hash_val = (hash_val << 2) | 0; break;
    //                     case 'T': hash_val = (hash_val << 2) | 1; break;
    //                     case 'C': hash_val = (hash_val << 2) | 2; break;
    //                     case 'G': hash_val = (hash_val << 2) | 3; break;
    //                     default: break;
    //                 }
    //             }
    //             if (hash_val < num_kmers) { // 防止越界
    //                 atomicAdd(&shared_counts[hash_val], 1);
    //             }
    //         }
    //     __syncthreads();

    //     for (int i = threadIdx.x; i < num_kmers; i += blockDim.x) {
    //         atomicAdd(&k_mer_Matrix[i], shared_counts[i]);
    //     }
    // }

    // __host__ void one_cluster::genClusters_GPU() {

    //     for (int i = 0 ;i < this->total_SeqList.size() ; i++) {
    //         this->hashMap[this->total_SeqList[i]] = false;
    //     }

    //     int seq_count = this->total_SeqList.size();

    //     while (seq_count >= 1) {
    //         int long_id = getLoneID();
    //         if (long_id >= 0) {
    //             this->cluster[long_id] = get_smi_seqs_GPU(long_id);
    //             seq_count -= this->cluster[long_id].size() + 1;
    //             this->hashMap[long_id] = true;
    //         }else {
    //             PLOGD << "==== 数组越界！====";
    //             exit(0);
    //         }
    //     }
    // }
    // std::vector<int> one_cluster::get_smi_seqs_GPU(int long_id) {
    //     std::vector<int> res;
    //     std::vector<int> kmerA = All_seqs[long_id].k_mer_list;
    //     for (int i = 0; i < total_SeqList.size(); i++) { // 从 0 开始会包含 long_id 自身
    //         if (!this->hashMap[this->total_SeqList[i]]) {
    //             std::vector<int> kmerB = All_seqs[this->total_SeqList[i]].k_mer_list;
    //             if (this->getDistanceCUDA(kmerA, kmerB) >= this->sim) {
    //                 res.push_back(this->total_SeqList[i]);
    //                 this->hashMap[this->total_SeqList[i]] = true;
    //             }
    //         }
    //     }
    //     return res;
    // }