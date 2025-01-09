//
// Created by aohy on 24-12-25.
//

#include "cluster.h"

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
}