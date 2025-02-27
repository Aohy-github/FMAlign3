//
// Created by aohy on 24-12-27.
//

#ifndef ALIGN_H
#define ALIGN_H

#include "top.h"
#include "bindings/cpp/WFAligner.hpp"


typedef struct Align_Seq
{
    std::string content;   
    // 保存所比对的中心序列的间隔
    std::vector<std::pair<int,int>> center_gap_list;
    // 保存自己的的间隔
    std::vector<std::pair<int,int>> self_gap_list;
    Align_Seq(std::string seq):content(seq) {};
}Align_Seq;


typedef struct Process_seq{
    int seq_id;
    std::string content;
    Process_seq(int id , std::string seq):seq_id(id),content(seq) {};
}Process_seq;


void mode_pick(std::vector<Process_seq>& process_seq_list);
std::vector<std::vector<std::pair<int,int>>> build_index_and_chain(std::vector<Process_seq> &process_seq_list , int k_mer);
void process_chain_table(std::vector<std::vector<std::pair<int,int>>>& chain_table);
std::vector<std::vector<Align_Seq>> get_split_chain_table(std::vector<std::vector<std::pair<int,int>>>& chain_table,std::vector<Process_seq>& process_seq_list , int k_mer);
std::vector<std::vector<std::pair<int, int>>> create_chain_(std::vector<std::vector<std::pair<int, int>>> chain_table , std::vector<Process_seq>& process_seq_list,int k_mer);
bool merge_split(std::vector<std::vector<Align_Seq>>& split_seq_list , std::vector<Process_seq>& process_seq_list);

int get_best_location_ssw(int main_location,int center_id , std::string& center_content , int querry_location , int seq_id , std::string& querry_content,int k_mer);

void process_cigar(std::string cigar,std::string& patter, std::vector<std::pair<int,int>>& gap_center_list);

void align_wfa_PSA(std::string & text , std::string & pattern, std::vector<std::pair<int,int>>& gap_center_list);

std::vector<std::pair<int, int>> get_lis(const std::vector<std::pair<int, int>>& data);

void cout_chain_table(std::vector<std::vector<std::pair<int,int>>>& chain_table);

void star_align(std::vector<std::vector<Align_Seq>>& split_seq_list);


void align_wfa_profile(std::string & text , std::string & pattern);

void process_cigar_two(std::string cigar, std::string& pattern , std::string& res_patter,std::string& text, std::string& res_text);

std::vector<std::pair<int,int>> merge(std::vector<std::pair<int,int>> gap_list);

void insert_gap(std::vector<std::pair<int,int>>& gap_list , std::string& pattern);

void filter_gap(std::vector<std::pair<int,int>>& patter1_gap_list, 
               std::vector<std::pair<int,int>>& all_merged_gaps, 
               std::vector<std::pair<int,int>>& center_gap_list);

#endif //ALIGN_H
