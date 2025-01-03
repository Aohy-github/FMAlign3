//
// Created by aohy on 24-12-24.
//

#ifndef TOP_H
#define TOP_H

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cuda_runtime_api.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <map>
// #include <thrust/device_vector.h>
// #include <thrust/transform_reduce.h>
// #include <thrust/functional.h>
#include <plog/Log.h> // 必须包含
#include <plog/Initializers/RollingFileInitializer.h>
#include <algorithm>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_pipeline.h>
#include <tbb/global_control.h>
#include <tbb/task_arena.h>
#include <tbb/concurrent_vector.h>
#include <tbb/tbb.h>
#include <sys/stat.h>
#include "bwt.h"
#include "ssw_cpp.h"
#include "ssw.h"





inline int K = 6;

typedef struct Sequences
{
    std::string name;
    std::string content;
    long long lens = 0;
    int k_mer_size = 0;
    long long A_A; // 内积
    std::vector<int> k_mer_list;
    // 总的k-mer可能性为 4^K
    Sequences(std::string name, std::string seq, int K , int lens):name(name), content(seq), k_mer_size(K), k_mer_list(1 << (2*K),0) , lens(lens) {};

}Seq;

inline std::vector<Seq> All_seqs;





#include "readfile.h"
#include "fcgr.h"
#include "cluster.h"
#include "center_align.h"
#include "fmindex.h"
#endif //TOP_H
