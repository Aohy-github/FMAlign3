//
// Created by aohy on 24-12-24.
//

#ifndef READFILE_H
#define READFILE_H
#include <top.h>

std::vector<Seq> get_sequence(const std::string& filename);
bool out_sequence(const std::string& filename);
#endif 
