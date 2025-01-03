//
// Created by aohy on 24-12-27.
//

#ifndef ALIGN_H
#define ALIGN_H

#include "top.h"
#include "bindings/cpp/WFAligner.hpp"

void center_align(std::vector<int> seq_id_list);
void build_index_and_chain(std::vector<int> seq_id_list , int k_mer);
void process_chain(std::vector<std::vector<std::pair<int,int>>>& chain_table , std::vector<int>& cover_table , int k_mer , std::vector<int> seq_id_list);




int get_best_location_ssw(int main_location,int center_id , int querry_location , int seq_id,int k_mer, double spand);
int process_cigar(std::string cigar, std::string& pattern , std::string& res_patter, std::vector<int>& res_text);
std::vector<int> align_wfa(std::string & text , std::string & pattern);
std::vector<std::pair<int, int>> get_lis(const std::vector<std::pair<int, int>>& data);
void cout_chain_table(std::vector<std::vector<std::pair<int,int>>>& chain_table , std::vector<int>& cover_table);
void process_more_seq_and_lens_less_10000(std::vector<int> seq_id_list);

#endif //ALIGN_H
