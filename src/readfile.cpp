

#include "readfile.h"



std::vector<Seq> get_sequence(const std::string& filename){
    int max_lens = 0, min_lens = INT32_MAX , mid_lens = 0;
    std::ifstream infile(filename);
    std::vector<Seq> seqs;
    if(!infile.is_open()) {
        printf("open file failed!");
        exit(0);
    }
    else {
        std::string line;
        std::string name;
        std::string content;
        while(getline(infile, line)){
            if(line[0] == '>') {
                if(content.size() > 0) {
                    max_lens = std::max(max_lens, (int)content.size());
                    min_lens = std::min(min_lens, (int)content.size());
                    mid_lens += content.size();
                    Seq seq(name, content, K , content.size());
                    seqs.push_back(seq);
                }
                name = line;
                content.clear();
            }else {
                content += line;
            }
        }
        if(content.size() > 0) {
            Seq seq(name, content, K,content.size());
            max_lens = std::max(max_lens, (int)content.size());
            min_lens = std::min(min_lens, (int)content.size());
            mid_lens += content.size();
            seqs.push_back(seq);
        }
        mid_lens = mid_lens / seqs.size();
        infile.close();
    }
    PLOGD << "max lens :" << max_lens;
    PLOGD << "min lens :" << min_lens;
    PLOGD << "mid lens :" << mid_lens;
    std::sort(seqs.begin(), seqs.end(),[&](const Seq&a , const Seq&b) {
        // std::cout << "lens :" << a.lens << std::endl;
        return a.lens > b.lens; // 从大到小
    });
    return seqs;
}


bool out_sequence(const std::string& filename){
    std::ofstream outfile(filename);
    if(!outfile.is_open()) {
        printf("open file failed");
        exit(0);
    }
    for(auto& seq : All_seqs) {
        outfile << seq.name << std::endl;
        outfile << seq.content << std::endl;
    }
    outfile.close();
    return true;
}