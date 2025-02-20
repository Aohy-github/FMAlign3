//
// Created by aohy on 24-12-27.
//

#include "center_align.h"
#include <unordered_set>

void mode_pick(std::vector<int> seq_id_list) {
    PLOGD << "first id : " << seq_id_list[0];
    int seq_lens = All_seqs[seq_id_list[0]].content.size();
    int k_mer = 39;
    PLOGD << "seq_lens : " << seq_lens;
    auto time_start = std::chrono::high_resolution_clock::now();
    //build_index_and_chain(seq_id_list , k_mer);
    if(seq_lens >= 10000 && seq_id_list.size() > 2){   
        PLOGD <<"处理多条长度大余10000的序列";
        /*
            1.构建fmindex，构建分割点
            3.整合分割点
            3.分割文件，并行合并
        */
        std::vector<std::vector<std::pair<int,int>>> chain_table = build_index_and_chain(seq_id_list , k_mer);   
        std::vector<std::vector<MID_Seq>> split_seq_list = get_split_chain_table(chain_table,seq_id_list, k_mer);
        
        star_align(seq_id_list , split_seq_list);
        merge_split(split_seq_list ,seq_id_list);
        std::cout << "==== FINISH BUILD INDEX AND CHAIN ====" << std::endl;
    }
    else if(seq_lens <= 10000 && seq_id_list.size() >= 2){ 
        PLOGD<< "处理多条长度小余10000的序列";

        /**
         * 1.进行中心序列的比对
         */
        std::vector<std::vector<MID_Seq>> not_split_seq_list;
        std::vector<MID_Seq> tmp;
        for(int i = 0; i < seq_id_list.size(); i++){
            tmp.push_back(MID_Seq(All_seqs[seq_id_list[i]].content));
        }
        not_split_seq_list.push_back(tmp);
        star_align(seq_id_list , not_split_seq_list);
        merge_split(not_split_seq_list,seq_id_list);
        // 输出中间结果
        // for(int i = 0; i < seq_id_list.size(); i++){
        //     PLOGD <<" , seq_name : " << All_seqs[seq_id_list[i]].name << ": " << All_seqs[seq_id_list[i]].content;
        // }
    }else if(seq_lens >= 10000 && seq_id_list.size() == 2){ 
        PLOGD << "处理两条序列 但很长的序列";
        std::vector<std::vector<std::pair<int,int>>> chain_table = build_index_and_chain(seq_id_list , k_mer);
        std::vector<std::vector<MID_Seq>> split_seq_list = get_split_chain_table(chain_table,seq_id_list, k_mer);
        
        star_align(seq_id_list , split_seq_list);
        merge_split(split_seq_list ,seq_id_list);   

    }else if(seq_lens <= 10000 && seq_id_list.size() == 2){ 
        PLOGD << "处理两条 但不是很长的序列";
        std::vector<std::vector<MID_Seq>> not_split_seq_list;

        for(int i = 0; i < seq_id_list.size(); i++){
            std::vector<MID_Seq> tmp;
            tmp.push_back(MID_Seq(All_seqs[seq_id_list[i]].content));
            not_split_seq_list.push_back(tmp);
        }

        star_align(seq_id_list , not_split_seq_list);
        merge_split(not_split_seq_list ,seq_id_list);
    }else{                                                 
        PLOGD << "序列数量太少，不比较";
    }
    auto time_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> time = time_end - time_start;
    std::cout << "==== FINISH MODE PICK ====" << "time : " << time.count() << " ms" << std::endl;
}



std::vector<std::vector<std::pair<int,int>>> build_index_and_chain(std::vector<int> seq_id_list , int k_mer){
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
    tmpfile[sizeof(tmpfile) - 1] = '\0';
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

    int step = k_mer / 2 > 15 ? 15 : k_mer / 2;
    int next_spand = step;
    bool bo = false;
    
    int table_size = (center_seq.size() - k_mer + 1) / step;
    std::vector<std::vector<std::pair<int,int>>> chain_table(seq_id_list.size(), std::vector<std::pair<int,int>>(table_size , {0,-1}));
    

    
    std::vector<unsigned int> locations_list;
    const int CENTER_SEQ_LENS = center_seq.size();
    for(uint64_t i = 0; i + k_mer < center_seq.size() - k_mer; i += next_spand){
        std::string patter_str = center_seq.substr(i, k_mer);
        
        
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
                       
        }
        if(seq_id_set.size() == 0){
            next_spand = k_mer + k_mer;
        }else{
            next_spand = k_mer / 2 > 15? 15 : k_mer / 2;
            for(auto it = seq_id_set.begin(); it != seq_id_set.end(); it++){
                chain_table[*it - 1][i / step] = {i / step , -1};
            }
        }
    }
    PLOGD << "====THIS TIME CREATE CHAIN FINISH ====";

    process_chain_table(chain_table , k_mer, seq_id_list);
    create_chain_(chain_table, seq_id_list,k_mer);
    std::cout << "==== FINISH CREATE CHAIN ====" << std::endl;
    return chain_table;
}
void process_chain_table(std::vector<std::vector<std::pair<int,int>>>& chain_table , int k_mer , std::vector<int> seq_id_list){
        
    PLOGD << "==== process_chain_table ====";
    cout_chain_table(chain_table);
    for(int i = 0 ; i < chain_table.size(); i++){
        for(int j = 0; j < chain_table[i].size(); j++){
        	PLOGD << "chain_table[" << i << "][" << j << "] : " << chain_table[i][j].first << " : " << chain_table[i][j].second;
        }
    }
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

    //std::cout << std::endl;
    //std::cout << std::endl;

    //cout_chain_table(chain_table);

    PLOGD << "==== FINISH process_chain_table ====";    
}
std::vector<std::vector<std::pair<int, int>>> create_chain_(std::vector<std::vector<std::pair<int, int>>> chain_table , std::vector<int> seq_id_list,int k_mer){
    
    std::cout << "==== create_chain_ ====" << std::endl;
    
    // 从上到下遍历，收集覆盖情况，收集两种情况，全覆盖和范围覆盖（0.9）
    std::vector<std::vector<std::pair<int,int>>> tmp_chain_table(chain_table.size());
    int cover_lens = seq_id_list.size();
    for(int i = 0; i < chain_table[0].size() ;i ++){ // 列遍历
        int tmp_cover_lens = 0;
        std::vector<std::pair<int,int>> tmp_complete_chain;
        
        if(chain_table[0][i].second == -1){
            continue;
        }
        for(int j = 0; j < chain_table.size() ; j++){ //行遍历
            if(chain_table[j][i].second != -1){
                tmp_cover_lens++;
                tmp_complete_chain.push_back(chain_table[j][i]); // 某一列的种子
            }else{
                int main_location = chain_table[0][i].second;
                int center_id = seq_id_list[0];
                int querry_location = main_location;
                int seq_id = seq_id_list[j];

                int best_location = get_best_location_ssw(main_location, center_id, querry_location, seq_id, k_mer);

                if(best_location != -1){
                    chain_table[j][i].second = best_location;
                    tmp_cover_lens++;
                }
                tmp_complete_chain.push_back(chain_table[j][i]);
            }
        }
        if(tmp_cover_lens >= cover_lens){
            for(int i = 0 ; i < tmp_complete_chain.size(); i++){
                tmp_chain_table[i].push_back(tmp_complete_chain[i]);
            }
        }
    }
    // 过滤 complete_chain_table 末尾的空值
    for(int j = tmp_chain_table[0].size() - 1; j >= 0; j-- ){
        if(tmp_chain_table[0][j].second == -1 || tmp_chain_table[0][j].first == 0){
            for(int i = 0; i < tmp_chain_table.size(); i++){
                tmp_chain_table[i].pop_back();
            }
        }else{
            break;
        }
    }
    process_chain_table(tmp_chain_table, k_mer, seq_id_list);
    cout_chain_table(tmp_chain_table);
    std::cout << "==== FINISH create_chain_ ====" << std::endl;
    return tmp_chain_table;
}
std::vector<std::vector<MID_Seq>> get_split_chain_table(std::vector<std::vector<std::pair<int,int>>>& chain_table,std::vector<int>& seq_id_list , int k_mer){
    PLOGD << "===== get_split_chain_table =====";
    std::vector<std::vector<std::pair<int,int>>> res_chain(seq_id_list.size());
    int pre = 0;
    int spand = 20000;
    for(int i = 0; i < chain_table[0].size() ; i++){ // 列遍历
        int location = chain_table[0][i].second;
        int distance = (location - pre);
        if(distance > spand){
            int count = 0;
            std::vector<std::pair<int,int>> tmp_chain;
            for(int j = 0; j < chain_table.size() ;j++){
                if(chain_table[j][i].second != -1){
                    count++;
                    tmp_chain.push_back(chain_table[j][i]);
                }
            }            
            if(count >= seq_id_list.size()){
                for(int j = 0; j < tmp_chain.size(); j++){
                    res_chain[j].push_back(tmp_chain[j]);
                }
                pre = location;
            }
        }
    }
    for(int i = 0; i < seq_id_list.size(); i++){
        res_chain[i].push_back({All_seqs[seq_id_list[i]].content.size(),All_seqs[seq_id_list[i]].content.size() - k_mer});
    }
    cout_chain_table(res_chain);

    std::vector<std::vector<MID_Seq>> split_seq_list;
    std::vector<int> pre_location(seq_id_list.size() , 0);
    std::vector<MID_Seq> tmp_split_list;
    // 第一个分割点
    for(int i = 0; i < res_chain.size() ; i ++){
        tmp_split_list.push_back(MID_Seq(All_seqs[seq_id_list[i]].content.substr(0,res_chain[i][0].second + k_mer)));
        pre_location[i] = res_chain[i][0].second;        
    }
    split_seq_list.push_back(tmp_split_list);
    //后续分割点
    for(int j = 1; j < res_chain[0].size(); j ++){
        tmp_split_list.clear();
        for(int i = 0; i < res_chain.size() ; i ++){
            tmp_split_list.push_back(MID_Seq(All_seqs[seq_id_list[i]].content.substr(pre_location[i] + k_mer,res_chain[i][j].second - pre_location[i])));
            pre_location[i] = res_chain[i][j].second;
        }
        split_seq_list.push_back(tmp_split_list);
    }
    // std::vector<int> seq_lens(seq_id_list.size() , 0);
    // for(int i = 0; i < split_seq_list.size(); i++){
    //     for(int j = 0; j < split_seq_list[i].size(); j++){
    //         seq_lens[i] += split_seq_list[i][j].content.size();
    //     }
    // }
    // for(int i = 0 ; i < seq_lens.size(); i++){
    //     std::cout << "seq_id : " << seq_id_list[i] << ", seq_len : " << seq_lens[i] << std::endl;
    // }
    for(int i = 0; i < seq_id_list.size(); i++){
        All_seqs[seq_id_list[i]].content.clear();
    }
    return split_seq_list;
}
void cout_chain_table(std::vector<std::vector<std::pair<int,int>>>& chain_table){
    int column_width = 7;
    for(int time = 0; time * 20 < chain_table[0].size(); time++)
    {
        for(int i = 0 ; i < chain_table.size() ; i ++){
            int count = 20;
            for(int j = time * 20; j < chain_table[i].size() && count ; j++){
                std::cout <<std::right
                          <<std::setw(column_width)
                          << chain_table[i][j].first 
                          << ":"<<std::setw(column_width)
                          << chain_table[i][j].second << ", ";
                count--;
            }
            std::cout << std::endl;

        }
        std::cout << std::endl;
    }
}

int get_best_location_ssw(int main_location,int center_id , int querry_location , int seq_id,int k_mer) {
    // 保证尾部对齐性
    // Declares a default Aligner
    // match_score, mismatch_penalty, gap_opening_penalty, gap_extending_penalty
    StripedSmithWaterman::Aligner aligner(2,2,2,1);
    // Declares a default filter
    StripedSmithWaterman::Filter filter;
    // Declares an alignment that stores the result

    StripedSmithWaterman::Alignment alignment;
    
    int sorce = (k_mer - 2) * 2; // 可接受错误的分数
    // 调用 ssw 计算是否存在相似k-mer
    //PLOGD << "main_location : " << main_location << ", center_id : " << center_id << ", querry_location : " << querry_location << ", seq_id : " << seq_id << ", k_mer : " << k_mer;
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

        //PLOGD << "querry_seq_kmer : " << querry_seq_kmer;
        //PLOGD << "ref_main_kmer : " << ref_main_kmer;

        //aligner.Align( ref_main_kmer.c_str(),querry_seq_kmer.c_str(), querry_seq_kmer.size(), filter, &alignment, maskLen);
        aligner.Align(querry_seq_kmer.c_str(), ref_main_kmer.c_str(), ref_main_kmer.size(), filter, &alignment, maskLen);
        std::string cigar = alignment.cigar_string;
        //PLOGD << "cigar_string : " << cigar;
        //PLOGD << "sw_score : " << alignment.sw_score;
        //PLOGD << "sorce : " << sorce;
        if(alignment.sw_score >= sorce ) {
            int op = 0;
            int num = 0;
            int location = sub_seq_querry_location + alignment.ref_begin;
            while(op < cigar.size()){
                if(isdigit(cigar[op])){
                    num = num * 10 + cigar[op] - '0';
                }else{
                    if(cigar[op] == 'S'){
                        location -= num;
                        break;
                    }else{
                        break;
                    }
                }
                op++;
            }
            if(cigar_int_to_op(cigar[0]) == 'S'){
                num = cigar_int_to_len(cigar[0]);
                location -= num;
            }
            if(location >= 0){
                //PLOGD << "match k-mer :" << All_seqs[seq_id].content.substr(location, k_mer);
                return location;
            }else{
                return -1;
            }
            
        }else{
            return -1;
        }
    }else{
        return -1;
    }
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

void align_wfa_PSA(std::string & text , std::string & pattern, std::vector<std::pair<int,int>>& center_gap_list){

    int text_size = text.size();
 
    wfa::WFAlignerGapAffine aligner(2,3,1,wfa::WFAligner::Alignment,wfa::WFAligner::MemoryHigh); 
    // wfa::WFAlignerGapLinear aligner(3,1,wfa::WFAligner::Alignment,wfa::WFAligner::MemoryHigh); 
    aligner.alignEnd2End(pattern,text);
    process_cigar(aligner.getAlignment(), pattern,center_gap_list);
    
}
void star_align(std::vector<int> seq_id_list , std::vector<std::vector<MID_Seq>>& split_seq_list){
    PLOGD << "process_more_seq_and_lens_less_10000 : seq_id_list";
    std::vector<std::pair<int,int>> All_gap_list;

    // for(int i = 0; i < split_seq_list.size(); i++){
    // 	std::vector<MID_Seq> tmp = split_seq_list[i];
    //     for(int j = 0; j < tmp.size(); j++){
    //     	PLOGD << tmp[j].content;
    //     }
    //     PLOGD << " ===================== ";
    // }
    // 设置全局线程数限制，例如 4 个线程
    PLOGD << "split_seq_list.size() : " << split_seq_list.size();
    for(int i = 0; i < split_seq_list.size(); i++){
        PLOGD << "i : " << i;
        std::vector<MID_Seq>& mid_seq = split_seq_list[i];
        std::string text_seq = mid_seq[0].content;
        PLOGD << "text_seq.size() : " << text_seq.size();
        for(int j = 1; j < mid_seq.size(); j++){
            PLOGD << "content.size() : " << mid_seq[j].content.size();
            std::vector<std::pair<int,int>> tmp_center_list;
            align_wfa_PSA(text_seq, mid_seq[j].content, tmp_center_list);
            mid_seq[0].center_gap_list.insert(mid_seq[0].center_gap_list.end(), tmp_center_list.begin(), tmp_center_list.end());
            mid_seq[j].center_gap_list.insert(mid_seq[j].center_gap_list.end(), tmp_center_list.begin(), tmp_center_list.end());
        }
        // 合并所有的中心序列gap ,合并所有交集
        std::vector<std::pair<int,int>> center_merged_gaps = merge(mid_seq[0].center_gap_list);
        //对中心序列进行插空
        insert_gap(center_merged_gaps, mid_seq[0].content);

        std::sort(center_merged_gaps.begin(), center_merged_gaps.end(), [](const std::pair<int,int>& a, const std::pair<int,int>& b){
            return a.first < b.first;
        });
        // 对其余序列进行插空
        for(int j = 1; j < mid_seq.size(); j++){
            filter_gap(mid_seq[j].self_gap_list,center_merged_gaps, mid_seq[j].center_gap_list);
            insert_gap(mid_seq[j].self_gap_list, mid_seq[j].content);
        }
    }
    // tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, 2);

    // tbb::parallel_for(tbb::blocked_range<int>(0, split_seq_list.size(), 5), [&split_seq_list](const tbb::blocked_range<int>& r) {
    //     for(int i = r.begin(); i < r.end(); i++){
    //         std::vector<MID_Seq>& mid_seq = split_seq_list[i];
    //         std::string text_seq = mid_seq[0].content;
    //         for(int j = 1; j < mid_seq.size(); j++){
    //             std::vector<std::pair<int,int>> tmp_center_list;
    //             align_wfa_PSA(text_seq, mid_seq[j].content, tmp_center_list);
    //             mid_seq[0].center_gap_list.insert(mid_seq[0].center_gap_list.end(), tmp_center_list.begin(), tmp_center_list.end());
    //             mid_seq[j].center_gap_list.insert(mid_seq[j].center_gap_list.end(), tmp_center_list.begin(), tmp_center_list.end());
    //         }
    //         // 合并所有的中心序列gap ,合并所有交集
    //         std::vector<std::pair<int,int>> center_merged_gaps = merge(mid_seq[0].center_gap_list);
    //         //对中心序列进行插空
    //         insert_gap(center_merged_gaps, mid_seq[0].content);

    //         std::sort(center_merged_gaps.begin(), center_merged_gaps.end(), [](const std::pair<int,int>& a, const std::pair<int,int>& b){
    //             return a.first < b.first;
    //         });
    //         // 对其余序列进行插空
    //         for(int j = 1; j < mid_seq.size(); j++){
    //             filter_gap(mid_seq[j].self_gap_list,center_merged_gaps, mid_seq[j].center_gap_list);
    //             insert_gap(mid_seq[j].self_gap_list, mid_seq[j].content);
    //         }
    //         // split_seq_list[i] = mid_seq;
    //     }
    // });
    // tbb::parallel_for_each(split_seq_list.begin(), split_seq_list.end(),[&](std::vector<MID_Seq>& mid_seq){
    //     int count = 0;
    //     std::string text_seq = mid_seq[0].content;
        
    //     for(int i = 1 ; i < mid_seq.size(); i++){
    //         std::cout << count ++ << std::endl;
    //         std::vector<std::pair<int,int>> tmp_center_list;
    //         align_wfa_PSA(text_seq, mid_seq[i].content, tmp_center_list);
    //         mid_seq[0].center_gap_list.insert(mid_seq[0].center_gap_list.end(), tmp_center_list.begin(), tmp_center_list.end());
    //         mid_seq[i].center_gap_list.insert(mid_seq[i].center_gap_list.end(), tmp_center_list.begin(), tmp_center_list.end());
    //     }
    //     // 合并所有的中心序列gap ,合并所有交集
    //     std::vector<std::pair<int,int>> center_merged_gaps = merge(mid_seq[0].center_gap_list);
    //     //对中心序列进行插空
    //     insert_gap(center_merged_gaps, mid_seq[0].content);

    //     std::sort(center_merged_gaps.begin(), center_merged_gaps.end(), [](const std::pair<int,int>& a, const std::pair<int,int>& b){
    //         return a.first < b.first;
    //     });
    //     // 对其余序列进行插空
    //     for(int i = 1; i < mid_seq.size(); i++){
    //         filter_gap(mid_seq[i].self_gap_list,center_merged_gaps, mid_seq[i].center_gap_list);
    //         insert_gap(mid_seq[i].self_gap_list, mid_seq[i].content);
    //     }
        
    // });

}

void process_cigar_two(std::string cigar, std::string& pattern , std::string& res_patter,std::string& text, std::string& res_text){
    int count_X = 0 , count_I = 0, count_M = 0 , count_D = 0;
    int op_pa = 0 , op_te = 0;
    for(int i = 0 ; i < cigar.size(); i++){
        char ch = cigar[i];
        switch(ch){
            case 'X' :{
                res_patter += pattern[op_pa++];
                res_text += text[op_te++];
                count_X++;
                break;
            }
            case 'I' :{
                res_patter += '-';
                count_I++;
                res_text += text[op_te++];
                break;
            }
            case 'M' :{
                res_patter += pattern[op_pa++];
                res_text += text[op_te++];
                count_M++;
                break;
            }
            case 'D' :{
                res_patter += pattern[op_pa++];
                res_text += '-';
                count_D++;
                break;
            }default:{
                std::cout << "other char : " << ch <<  std::endl;
                break;
            }
        }
    }
    
}

void align_wfa_profile(std::string & text , std::string & pattern){

    int text_size = text.size();
    
    std::string res_patter = "";
    std::string res_text = "";
    
    PLOGD << "text : " << text;
    PLOGD << "pattern : " << pattern;
    
    wfa::WFAlignerGapAffine aligner(4,6,2,wfa::WFAligner::Alignment,wfa::WFAligner::MemoryHigh); // 进行对齐需要多少消耗
    
    aligner.alignEnd2End(pattern,text);
    process_cigar_two(aligner.getAlignment(), pattern,res_patter,text,res_text);

    pattern = res_patter;
    text = res_text;
}

std::string get_common_seq(std::vector<int> seq_id_list){
    std::vector<std::vector<int>> count_alph;
    count_alph.resize(All_seqs[seq_id_list[0]].content.size(), std::vector<int>(5,0));
    for(int i = 0; i < seq_id_list.size(); i++){
        for(int j = 0; j < All_seqs[seq_id_list[i]].content.size(); j++){
            char ch = All_seqs[seq_id_list[i]].content[j];
            switch(ch){
                case 'A':{
                    count_alph[j][0]++;
                    break;
                }
                case 'T':{
                    count_alph[j][1]++;
                    break;
                }
                case 'C':{
                    count_alph[j][2]++;
                    break;
                }
                case 'G':{
                    count_alph[j][3]++;
                    break;
                }
                case 'a':{
                    count_alph[j][0]++;
                    break;
                }
                case 't':{
                    count_alph[j][1]++;
                    break;
                }
                case 'c':{
                    count_alph[j][2]++;
                    break;
                }
                case 'g':{
                    count_alph[j][3]++;
                    break;
                }
                case '-':{
                    count_alph[j][4]++;
                    break;
                }
            }
        }
    }
    std::string final_seq = "";
    for(int i = 0; i < count_alph.size(); i++){
        int max = 0;
        int index = 0;
        int count = seq_id_list.size();
        for(int j = 0; j < count_alph[i].size(); j++){
            if(count_alph[i][j] > max){
                max = count_alph[i][j];
                index = j;
            }
        }
        if(max >= (int)(count * 0.9)){
            switch(index){
                case 0:{
                    final_seq += "A";
                    break;
                }
                case 1:{
                    final_seq += "T";
                    break;
                }
                case 2:{
                    final_seq += "C";
                    break;
                }
                case 3:{
                    final_seq += "G";
                    break;
                }
            }

        }else{
            final_seq += "N";
        }
    }
    int pre = -1;
    for(int i = 0; i < final_seq.size(); i++){
        if(final_seq[i] == 'N'){
            if(pre == -1){
                pre = i;
                continue;
            }else{
                if(i - pre < 15){
                    for(int j = pre; j < i; j++){
                        final_seq[j] = 'N';
                    }
                }
            }
        }
    }
    return final_seq;
}

std::vector<std::pair<int,int>> merge(std::vector<std::pair<int,int>> gap_list){ // 合并间隙，可能存在重复的间隙，保留同位置最大间隙
    if (gap_list.size() == 0) {
        return {};
    }
    //std::cout << "=== MERGE ==="<<std::endl;
    std::vector<std::pair<int,int>> merged;
    std::unordered_map<int,int> final_gap_list;
    for(auto gap : gap_list){
        int start = gap.first;
        int lens = gap.second;
        if(final_gap_list.find(start) == final_gap_list.end()){
            final_gap_list[start] = lens;
        }else{
            final_gap_list[start] = std::max(final_gap_list[start], lens);
        }
    }
    
    
    for(auto gap : final_gap_list){
        //std::cout << gap.first << " " << gap.second << std::endl;
        merged.emplace_back(std::pair<int,int>{gap.first, gap.second}); 
    }
    return merged;
}
void process_cigar(std::string cigar,std::string& patter, std::vector<std::pair<int,int>>& gap_center_list){
    int count_center = 0;
    int op_pa = 0 , op_te = 0;
    std::string res_patter = "";
    for(char ch : cigar){
        // D 对text插空， I 对pattern插空， M 匹配， X 不匹配
        if(ch == 'D'){
          count_center++;
          res_patter += patter[op_pa++];
        }
        else if(ch == 'I'){
          res_patter += '-';
          op_te++;
        }
        else if(ch =='M' || ch == 'X'){ // 'M' 或 'X'
            
            // 处理连续的 'D' 操作
            if(count_center > 0){
                gap_center_list.emplace_back(std::pair<int,int>{op_te, count_center});
                count_center = 0;
            }
            // 处理当前的 'M' 或 'X' 操作
            op_te++;
            res_patter += patter[op_pa++];
        }else{
            std::cout << "other char : " << ch << std::endl;
            exit(0);
        }
    }

    // 处理尾部剩余的 'D' 操作

    if(count_center > 0){
        gap_center_list.emplace_back(std::pair<int,int>{op_te, count_center});
    }
    // std::cout << "gap_center_list.size():" << gap_center_list.size() << std::endl;
    // for(auto gap: gap_center_list){
    //     std::cout << gap.first << " " << gap.second << std::endl;
    // }
    patter = res_patter;
}

void insert_gap(std::vector<std::pair<int,int>>& gap_list , std::string& pattern){
    // 按照从高到低的顺序插入间隙
    sort(gap_list.begin(), gap_list.end(), [](const std::pair<int,int>& a, const std::pair<int,int>& b) {
        return a.first > b.first;
    });
    
    for(auto &gap : gap_list){
        int start = gap.first;
        int lens = gap.second;
        // 在 start 位置插入 (end - start) 个 '-'
        //std::cout << "start : " << start << " , end : " << lens << std::endl;
        pattern.insert(start, lens, '-');
    }
}


void filter_gap(std::vector<std::pair<int,int>>& patter1_gap_list, 
               std::vector<std::pair<int,int>>& all_merged_gaps, 
               std::vector<std::pair<int,int>>& center_gap_list) {
    // std::cout<< "=== filter_gap ====" <<std::endl;
    // std::cout<< "all_merged_gaps.size():" << all_merged_gaps.size() << std::endl;
    // for(int i = 0; i < all_merged_gaps.size();i++){
    //   std::cout << all_merged_gaps[i].first << " " << all_merged_gaps[i].second << std::endl;
    // }
    
    // std::cout<< "center_gap_list.size() :" << center_gap_list.size() << std::endl;
    sort(center_gap_list.begin(), center_gap_list.end(), [](const std::pair<int,int>& a, const std::pair<int,int>& b) {
        return a.first < b.first;
    });
    // for(int i = 0; i < center_gap_list.size();i++){
    //   std::cout << center_gap_list[i].first << " " << center_gap_list[i].second << std::endl;
    // }

    std::vector<std::pair<int,int>> tmp;

    size_t i = 0; // 指向 all_merged_gaps
    size_t j = 0; // 指向 center_gap_list
  

  /**
   * 
   * 1. all.first < self.first
   *      需要添加到tmp中
   * 2.all.first == self.first, 同时all.second >= self.second
   * 3.all.first > self.first
   * 
   * 
   * 
   */
    int offset = 0;
    while(i < all_merged_gaps.size() && j < center_gap_list.size()){
        std::pair<int,int> all = all_merged_gaps[i];
        std::pair<int,int> self= center_gap_list[j];

        if(all.first < self.first){
            tmp.push_back({all.first + offset, all.second});
            i++;
        }else if(all.first == self.first){
          int len = all.second - self.second;
          if(len > 0){
            tmp.push_back(std::pair<int,int>{all.first + offset, len});
          }
          offset += self.second;
          i++ , j++;
        }else{ // all.first > self.first
            tmp.push_back(std::pair<int,int>{all.first + offset, all.second});
            i++;
        }
    }
    while(i < all_merged_gaps.size()){
      tmp.push_back({all_merged_gaps[i].first + offset, all_merged_gaps[i].second});
      i++;
    }
    for(int i = 0; i < tmp.size();i++){
      //std::cout << tmp[i].first << " " << tmp[i].second << endl;
      patter1_gap_list.push_back(tmp[i]);
    }
}

bool merge_split(std::vector<std::vector<MID_Seq>>& split_seq_list , std::vector<int> seq_id_list){
    PLOGD << "===== merge_and_split =====";
    PLOGD << "seq_id_list.size() : " << seq_id_list.size();
    PLOGD << "split_seq_list.size() : " << split_seq_list.size();
    PLOGD << "split_seq_list[0].size() : " << split_seq_list[0].size();
    PLOGD << "All_seqs.size() : " << All_seqs.size();
    
    for(int i = 0 ; i < seq_id_list.size(); i++){
        All_seqs[seq_id_list[i]].content.clear();
    }
    for(int i = 0; i < split_seq_list.size(); i++){
    	std::vector<MID_Seq> tmp = split_seq_list[i];
        for(int j = 0; j < tmp.size(); j++){
            PLOGD << "id :" << seq_id_list[j] << " , content : " << tmp[j].content;
            All_seqs[seq_id_list[j]].content += tmp[j].content;
        }
        //PLOGD << " ===================== ";
    }
    PLOGD << "All_seqs.size() : " << All_seqs.size();
    return true;
}

