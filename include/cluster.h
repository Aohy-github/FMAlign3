//
// Created by aohy on 24-12-25.
//

#ifndef CLUSTER_H
#define CLUSTER_H

#include "top.h"
namespace Cluster {
    class Clu {
        public:
            Clu() {};
            ~Clu() {};
            // 总的聚类结果
            std::map<int,std::vector<int>> TOP_clusters;
            // 相对中心序列的距离
            std::map<int,double> to_center_distance;
            std::map<int,std::vector<int>> TOP_clusters_sort;
        int selectLongest();
        
        void get_cluster_two();
        void get_align_sort();
        std::vector<int> cluster_for_lens(int longest_id);
    };

    std::string get_common_seq(std::vector<std::string> content_list);
    void refresh_seq_content(std::vector<int>& seq_id_list , std::string seq1);
} // Cluster

#endif //CLUSTER_H
