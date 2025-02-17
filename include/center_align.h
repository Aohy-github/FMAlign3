//
// Created by aohy on 24-12-27.
//

#ifndef ALIGN_H
#define ALIGN_H

#include "top.h"
#include "bindings/cpp/WFAligner.hpp"

void mode_pick(std::vector<int> seq_id_list);
std::vector<std::vector<std::pair<int,int>>> build_index_and_chain(std::vector<int> seq_id_list , int k_mer);
void process_chain_table(std::vector<std::vector<std::pair<int,int>>>& chain_table , int k_mer , std::vector<int> seq_id_list);
std::vector<std::vector<MID_Seq>> get_split_chain_table(std::vector<std::vector<std::pair<int,int>>>& chain_table,std::vector<int>& seq_id_list , int k_mer);
std::vector<std::vector<std::pair<int, int>>> create_chain_(std::vector<std::vector<std::pair<int, int>>> chain_table , std::vector<int> seq_id_list,int k_mer);
bool merge_and_split(std::vector<std::vector<MID_Seq>>& split_seq_list , std::vector<int> seq_id_list);

int get_best_location_ssw(int main_location,int center_id , int querry_location , int seq_id,int k_mer);

void process_cigar(std::string cigar,std::string& patter, std::vector<std::pair<int,int>>& gap_center_list);

void align_wfa_PSA(std::string & text , std::string & pattern, std::vector<std::pair<int,int>>& gap_center_list);

std::vector<std::pair<int, int>> get_lis(const std::vector<std::pair<int, int>>& data);

void cout_chain_table(std::vector<std::vector<std::pair<int,int>>>& chain_table);

void star_align(std::vector<int> seq_id_list , std::vector<std::vector<MID_Seq>>& split_seq_list);

void process_seq_clusters_(std::map<int,std::vector<int>> clusters);

std::string get_common_seq(std::vector<int> seq_id_list);

void align_wfa_profile(std::string & text , std::string & pattern);

void process_cigar_two(std::string cigar, std::string& pattern , std::string& res_patter,std::string& text, std::string& res_text);

std::vector<std::pair<int,int>> merge(std::vector<std::pair<int,int>> gap_list);

void insert_gap(std::vector<std::pair<int,int>>& gap_list , std::string& pattern);

void filter_gap(std::vector<std::pair<int,int>>& patter1_gap_list, 
               std::vector<std::pair<int,int>>& all_merged_gaps, 
               std::vector<std::pair<int,int>>& center_gap_list);

#endif //ALIGN_H
