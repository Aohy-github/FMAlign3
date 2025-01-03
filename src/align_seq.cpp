//
// Created by aohy on 24-12-27.
//

#include "center_align.h"
#include <unordered_set>

void center_align(std::vector<int> seq_id_list) {
    PLOGD << "first id : " << seq_id_list[0];
    int seq_lens = All_seqs[seq_id_list[0]].content.size();
    int k_mer = 6;
    PLOGD << "seq_lens : " << seq_lens;
    build_index_and_chain(seq_id_list , k_mer);
    if(seq_lens > 10000 && seq_id_list.size() > 2){   
        PLOGD <<"处理多条长度大余10000的序列";
        /*
            1.构建fmindex，构建分割点
            3.整合分割点
            3.分割文件，并行合并
        */
        build_index_and_chain(seq_id_list , k_mer);   

    }
    else if(seq_lens <= 10000 && seq_id_list.size() >= 2){ 
        PLOGD<< "处理多条长度小余10000的序列";

        /**
         * 1.进行中心序列的比对
         * 
         */
        process_more_seq_and_lens_less_10000(seq_id_list);
        
    }else if(seq_lens >= 10000 && seq_id_list.size() == 2){ 
        PLOGD << "处理两条序列 但很长的序列";
        build_index_and_chain(seq_id_list , k_mer);   



    }else if(seq_lens <= 10000 && seq_id_list.size() == 2){ 
        PLOGD << "处理两条 但不是很长的序列";




    }else{                                                 
        PLOGD << "序列数量太少，不比较";



    }
    
}



void build_index_and_chain(std::vector<int> seq_id_list , int k_mer){
    PLOGD << "==== build_index_and_chain ====";
    PLOGD << "seq_id_list.size() : " << seq_id_list.size();
    std::string file_path = "./tmp_file/";
    std::string file_name = "tmp_" + std::to_string(seq_id_list[0]) + '_' + std::to_string(seq_id_list.size()) + ".fasta";
    std::string final_name = file_path + file_name;
    std::ofstream outfile(final_name);
    if(!outfile.is_open()){
        char buffer[1024];
        PLOGD << "file_name : " << final_name;
        if (getcwd(buffer, sizeof(buffer)) != nullptr) {
            PLOGD <<"当前工作目录:" << buffer;
        }
        PLOGD << "open file failed";
        exit(0);
    }
    unsigned long long total_seq_size = 0;
    std::string seqs_content = "";
    std::vector<int> seq_list;
    seq_list.resize(seq_id_list.size() + 1); // 首位存放0
    seq_list[0] = 0;


    for(int i = 0; i < seq_id_list.size(); i++){
        seqs_content += All_seqs[seq_id_list[i]].content;
        total_seq_size += All_seqs[seq_id_list[i]].content.size();
        seq_list[i + 1] = total_seq_size;
    }
    
    char* total_seq = const_cast<char*>(seqs_content.data());
    
    char tmpfile[150];
    strncpy(tmpfile, final_name.c_str(), sizeof(tmpfile) - 1);
    indenpendent_creadte_index(total_seq_size , &total_seq,8,tmpfile);
    
    load_index(tmpfile);

    /**
     * @brief 一张二维表
     */
    //TATCACC
    std::string center_seq = All_seqs[seq_id_list[0]].content;
    unsigned int top = 0, bot = 0;
	unsigned int pre_top = 0, pre_bot = 0;
    unsigned int tmp_SA_length = 0;
    
    const char* str = center_seq.c_str();
    std::unordered_map<std::string,std::vector<unsigned int>> tmp_kmer;
    std::unordered_map<std::string,std::vector<bool>> tmp_kmer_id_flag;

    int step = k_mer / 2;
    
    int table_size = (center_seq.size() - k_mer + 1) / step;
    std::vector<std::vector<std::pair<int,int>>> chain_table(seq_id_list.size(), std::vector<std::pair<int,int>>(table_size , {0,-1}));
    std::vector<int> cover_table(table_size, 0);

    
    std::vector<unsigned int> locations_list;
    
    for(uint64_t i = 0; i < center_seq.size() - k_mer; i+= step){
        std::string patter_str = center_seq.substr(i, k_mer);
        // PLOGD <<"";
        // PLOGD <<"";
        PLOGD << "===== NEXT =====";
        
        
        unsigned int location_count = 0;
        unsigned int tmp_num_locations = 0;

        if(tmp_kmer.count(patter_str)) {
            location_count = tmp_kmer[patter_str].size();
            locations_list = tmp_kmer[patter_str];
            
        }
        else{
            char *patter = const_cast<char*>(patter_str.c_str());
            
            location_count = count(patter, k_mer, &top, &bot, &pre_top, &pre_bot);  
            if(location_count == 0) {
                continue;
            }

            locations_list.resize(location_count , 0);
            locate(patter, top, bot,pre_top, pre_bot, locations_list.data(), k_mer, &tmp_SA_length);


            tmp_kmer[patter] = locations_list;
        }


        int location_main = i;
        // PLOGD << "location_count : " << location_count;
        // PLOGD << " patter : " << patter_str;


        // 把主序列位置存储写如table中
        chain_table[0][i / step] = {i / step , i};
        int cover_count = 1;
        std::unordered_set<int> seq_id_set;
        for(int j = 0; j < seq_id_list.size(); j++){
            seq_id_set.insert(j + 1);
        }
        for(int j = 0; j < location_count; j++) {
            //PLOGD << "locations_list[j] : " << locations_list[j];   
            
            auto it = std::upper_bound(seq_list.begin(), seq_list.end(), locations_list[j]); // 第一个大余locations_list[j]的位置
            if(it != seq_list.end()){
                auto id = std::distance(seq_list.begin(), it);  // 返回的时第一个大余locations_list[j]的位置, 初始位置为0，

                if(seq_id_list[id - 1] == seq_id_list[0]){
                    // PLOGD << " *** main *** ";
                    // PLOGD << "id : " << seq_id_list[id - 1];
                    // PLOGD << " , k_mer : " << All_seqs[seq_id_list[0]].content.substr(location_main, k_mer);
                    // PLOGD << "location_main : " << location_main;
                    seq_id_set.erase(id);
                    continue;
                }else if(locations_list[j] + k_mer <= seq_list[id]){
                    int location_seq = locations_list[j] - seq_list[id - 1];
                    // PLOGD << " ++++ seq ++++ ";
                    // PLOGD << "id : " << seq_id_list[id - 1]; // id 表示从
                                              
                    // PLOGD << "location_seq : " << location_seq << " ,k-mer : " << All_seqs[seq_id_list[id - 1]].content.substr(location_seq, k_mer);
                    if(chain_table[id - 1][i / step].second == -1){
                        
                        chain_table[id - 1][i / step] = {i / step , location_seq};
                        cover_count++;
                        seq_id_set.erase(id);
                    }else if(chain_table[id - 1][i / step].second != -1){
                        // 保留相对位置最小的种子
                        int dif_1 = i - chain_table[id - 1][i / step].second;
                        dif_1 = abs(dif_1);
                        int dif_2 = i - location_seq;
                        dif_2 = abs(dif_2);

                        if(dif_1 > dif_2){
                            chain_table[id - 1][i / step] = {i / step , location_seq};
                        }

                    }
                }
            }else{
                PLOGD << "not found";
            }
                       
            cover_table[i / step] = cover_count;
        }
        for(auto it = seq_id_set.begin(); it != seq_id_set.end(); it++){
            chain_table[*it - 1][i / step] = {i / step , -1};
        }
        
    }
    PLOGD << "====THIS TIME CREATE CHAIN FINISH ====";

    process_chain(chain_table , cover_table, k_mer, seq_id_list);
    

}
void process_chain(std::vector<std::vector<std::pair<int,int>>>& chain_table , std::vector<int>& cover_table , int k_mer , std::vector<int> seq_id_list){
        
    cout_chain_table(chain_table, cover_table);
    

    // 处理交叉
    tbb::parallel_for_each(chain_table.begin(), chain_table.end(), [&](std::vector<std::pair<int,int>>& one_chain){
        std::vector<std::pair<int,int>> tmp_chain = get_lis(one_chain);
        int index = 0;
        for(int i = 0 ;i < one_chain.size();i++){
            if(one_chain[i].first == tmp_chain[index].first){
                one_chain[i] = tmp_chain[index];
                index++;
            }else{
                one_chain[i].second = -1;
            }
        }
        
    });

    std::cout << std::endl;
    std::cout << std::endl;

    cout_chain_table(chain_table, cover_table);
}

std::vector<std::pair<int, int>> get_lis(const std::vector<std::pair<int, int>>& data) {
    struct NodeWithIndex {
        int parent;         // 父节点的索引
        int value;          // 当前节点的值
        int original_index; // 在原始数据中的索引
    };
    // 第一步：过滤掉第二个元素为 -1 的对
    std::vector<std::pair<int, int>> filtered_data;
    for(const auto& p : data){
        if(p.second != -1){
            filtered_data.push_back(p);
        }
    }

    // 如果过滤后没有有效数据，返回空
    if(filtered_data.empty()){
        return {};
    }

    int n = filtered_data.size();
    std::vector<NodeWithIndex> nodes(n, NodeWithIndex{ -1, 0, 0 });

    // tails[i] 将存储长度为 i+1 的所有递增子序列中最后一个元素的最小值
    std::vector<int> tails; 
    // indices[i] 将存储长度为 i+1 的递增子序列中最后一个元素在 filtered_data 中的索引
    std::vector<int> indices; 
    // predecessors[i] 存储构成 LIS 时每个元素的前驱元素的索引
    std::vector<int> predecessors_filtered(n, -1); 

    for(int i = 0; i < n; ++i){
        // 二分查找当前值应插入的位置
        int pos = lower_bound(tails.begin(), tails.end(), filtered_data[i].second) - tails.begin();

        if(pos == tails.size()){
            tails.push_back(filtered_data[i].second);
            indices.push_back(i);
        }
        else{
            tails[pos] = filtered_data[i].second;
            indices[pos] = i;
        }

        if(pos != 0){
            predecessors_filtered[i] = indices[pos-1];
        }
    }

    // 重构 LIS
    std::vector<std::pair<int, int>> lis;
    int k = indices.empty() ? -1 : indices.back();
    while(k >=0){
        lis.emplace_back(filtered_data[k]);
        k = predecessors_filtered[k];
    }
    reverse(lis.begin(), lis.end());

    return lis;
}


int get_best_location_ssw(int main_location,int center_id , int querry_location , int seq_id,int k_mer, double spand) {
    
    // Declares a default Aligner
    // match_score, mismatch_penalty, gap_opening_penalty, gap_extending_penalty
    StripedSmithWaterman::Aligner aligner(2,1,1,1);
    // Declares a default filter
    StripedSmithWaterman::Filter filter;
    // Declares an alignment that stores the result

    StripedSmithWaterman::Alignment alignment;
    
    int sorce = (k_mer - 2) * 2; // 可接受错误的分数
    // 调用 ssw 计算是否存在相似k-mer
    PLOGD << "main_location : " << main_location << ", center_id : " << center_id << ", querry_location : " << querry_location << ", seq_id : " << seq_id << ", k_mer : " << k_mer << ", spand : " << spand ;
    std::string ref_main_kmer = All_seqs[center_id].content.substr(main_location, k_mer);

    int sub_seq_querry_location = querry_location - k_mer>= 0 ? querry_location - k_mer: 0;
    
    int sub_seq_querry_length = k_mer * 3;
    
    if(sub_seq_querry_location >= 0 && sub_seq_querry_location < All_seqs[seq_id].content.size()){
        int32_t maskLen = sub_seq_querry_length / 2;
        maskLen = maskLen < 15 ? 15 : maskLen;
        
        std::string querry_seq_kmer = All_seqs[seq_id].content.substr(sub_seq_querry_location, sub_seq_querry_length);
        std::string tmp = querry_seq_kmer;
        querry_seq_kmer = ref_main_kmer;
        ref_main_kmer = tmp;

        PLOGD << "querry_seq_kmer : " << querry_seq_kmer;
        PLOGD << "ref_main_kmer : " << ref_main_kmer;

        //aligner.Align( ref_main_kmer.c_str(),querry_seq_kmer.c_str(), querry_seq_kmer.size(), filter, &alignment, maskLen);
        aligner.Align(querry_seq_kmer.c_str(), ref_main_kmer.c_str(), ref_main_kmer.size(), filter, &alignment, maskLen);
        
        PLOGD << "cigar_string : " << alignment.cigar_string;
        PLOGD << "sw_score : " << alignment.sw_score;
        PLOGD << "sorce : " << sorce;
        if(alignment.sw_score >= sorce ) { 
            return sub_seq_querry_location + alignment.ref_begin;
        }else{
            return -1;
        }
    }else{
        return -1;
    }
}

int process_cigar(std::string cigar, std::string& pattern , std::string& res_patter, std::vector<int>& res_text){
    int count_X = 0 , count_I = 0, count_M = 0 , count_D = 0;
    int op_pa = 0 , op_te = 0;
    for(int i = 0 ; i < cigar.size(); i++){
        char ch = cigar[i];
        switch(ch){
            case 'X' :{
                res_patter += pattern[op_pa++];
                res_text.push_back(0);
                count_X++;
                break;
            }
            case 'I' :{
                res_patter += '-';
                count_I++;
                res_text.push_back(0);
                break;
            }
            case 'M' :{
                res_patter += pattern[op_pa++];
                res_text.push_back(0);
                count_M++;
                break;
            }
            case 'D' :{
                res_patter += pattern[op_pa++];
                res_text.push_back(1);
                count_D++;
                break;
            }default:{
                std::cout << "other char : " << ch <<  std::endl;
                break;
            }
        }
    }
    return count_X + count_I + count_D;
}



std::vector<int> align_wfa(std::string & text , std::string & pattern){

    int text_size = text.size();
    
    std::string res_patter = "";
    std::vector<int> res_text_list;
    PLOGD << "text : " << text;
    PLOGD << "pattern : " << pattern;
    
    wfa::WFAlignerGapAffine aligner(4,6,2,wfa::WFAligner::Alignment,wfa::WFAligner::MemoryHigh); // 进行对齐需要多少消耗
    
    aligner.alignEnd2End(pattern,text);
    process_cigar(aligner.getAlignment(), pattern,res_patter,res_text_list);

    pattern = res_patter;
    return res_text_list;
}
void cout_chain_table(std::vector<std::vector<std::pair<int,int>>>& chain_table , std::vector<int>& cover_table){
    for(int time = 0; time * 20 < cover_table.size(); time++)
    {
        for(int i = 0 ; i < chain_table.size() ; i ++){
            int count = 20;
            for(int j = time * 20; j < chain_table[i].size() && count ; j++){
                std::cout <<std::setw(3)<< chain_table[i][j].first << ":"<<std::setw(3)<< chain_table[i][j].second << ", ";
                count--;
            }
            std::cout << std::endl;

        }
        std::cout << std::endl;
    }
}


void process_more_seq_and_lens_less_10000(std::vector<int> seq_id_list){
    PLOGD << "process_more_seq_and_lens_less_10000";
    
    tbb::enumerable_thread_specific<int> local_max([]() -> int { return 0; });

    tbb::enumerable_thread_specific<std::vector<int>> total_center_list;
    tbb::parallel_for_each(seq_id_list.begin() , seq_id_list.end(),[&](const int id){
        if(id == seq_id_list[0]){
            return;
        }
        std::string center_seq = All_seqs[seq_id_list[0]].content;
        std::string querry_seq = All_seqs[id].content;
        std::vector<int>& tmp_center_list = total_center_list.local();

        tmp_center_list = align_wfa(center_seq, querry_seq);
        int current_size = tmp_center_list.size();
        if(current_size > local_max.local()){
            local_max.local() = current_size;
        }
        All_seqs[id].content = querry_seq;

    });
    int longest_tmp_center_list = local_max.combine([](const int a, const int b) -> int {return std::max(a, b);});

    std::vector<int> final_center_list(longest_tmp_center_list,0);
    for(auto& center_list : total_center_list){
        for(int i = 0; i < center_list.size(); i++){
            final_center_list[i] |= center_list[i];
        }
    }
    std::string final_center_seq = "";
    int op = 0;
    for(int i = 0; i < final_center_list.size(); i++){
        final_center_seq += final_center_list[i] ? '-' : All_seqs[seq_id_list[0]].content[op++];
    }
    All_seqs[seq_id_list[0]].content = final_center_seq;
    
}