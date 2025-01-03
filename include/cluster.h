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
        int selectLongest();
        //void get_cluster_one();
        void get_cluster_two();
        void get_cluster_three();
        std::vector<int> cluster_for_lens(int longest_id);
    };



} // Cluster

#endif //CLUSTER_H
