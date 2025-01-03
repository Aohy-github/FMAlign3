//
// Created by aohy on 24-12-27.
//

#ifndef FMINDEX_H
#define FMINDEX_H

#include "top.h"
#include "bwt.h"
// 序列是被从小到大排序了的
namespace FM {
// (使用FM-index构建公共种子，分割序列，并行比对，合并结果)

// 构建FM-index
        void buildFMIndex(std::string& text);

// 查找公共种子
        void findCommonSeeds();

// 分割序列
        void splitSequences();

// 并行比对
        void parallelAlignment();

// 合并结果
        void mergeResults();
class fmindex {
    public:
        // 存放所有序列的内容，根据传入序列集合的相对位置
        std::string all_seq_content;
        // 存放累加所有序列的长度，根据传入序列集合的相对位置
        std::vector<std::pair<int,int>> seq_lens_list;

        fmindex(std::map<int,std::vector<int>> cluster){
            // 根据传入的聚类集合，将所有序列的内容拼接起来
            int sum_lens = 0;
            for(auto cluster_id : cluster){
                for(auto seq_id : cluster_id.second){
                    all_seq_content += All_seqs[seq_id].content;
                    sum_lens += All_seqs[seq_id].lens;

                    seq_lens_list.push_back(std::make_pair(seq_id,sum_lens));
                }
            }

        };

        std::map<int,std::string> getAlignres(){
            std::map<int,std::string> align_res;

            FM::buildFMIndex(all_seq_content);



            return align_res;
        }

};

} // FM

#endif //FMINDEX_H
