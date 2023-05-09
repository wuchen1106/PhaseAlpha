#include <vector>
#include <map>
#include <iostream>

#include "TChain.h"

int main(int argc, char ** argv){
    TChain * ichain = new TChain("t","t");
    ichain->Add(argv[1]);

    std::vector<double> * start[4] = {0};
    std::vector<double> * width[4] = {0};
    std::vector<double> * height[4] = {0};

    for (int ch = 0; ch<4; ch++){
        ichain->SetBranchAddress(Form("start%c",'A'+ch),&start[ch]);
        ichain->SetBranchAddress(Form("width%c",'A'+ch),&width[ch]);
        ichain->SetBranchAddress(Form("height%c",'A'+ch),&height[ch]);
    }


    std::map<int,int> counter[4];
    std::map<int,double> totalWidth[4];
    Long64_t nEntries = ichain->GetEntries();
    for (int e = 0; e<nEntries; e++){
        ichain->GetEntry(e);
        for (int ch = 0; ch<4; ch++){
            if (!height[ch]) continue;
            double prev_time = -1;
            for (int i = 0; i<height[ch]->size(); i++){
                if (start[ch]->at(i)<2400) continue;
                if (height[ch]->at(i)>30){
                    double cur_time = start[ch]->at(i);
                    //if (prev_time>0&&fabs(cur_time-prev_time)<14){
                    //    std::cout<<"Event "<<e<<" ch"<<'A'+ch<<" hit "<<i<<" is too close to the previous one! dT = "<<cur_time<<" - "<<prev_time<<" = "<<cur_time-prev_time<<std::endl;
                    //}
                    prev_time = cur_time;
                }
                for (int thr = 10; thr<200; thr+=10){
                    if (height[ch]->at(i)>thr){
                        counter[ch][thr]++;
                        totalWidth[ch][thr]+=width[ch]->at(i);
                    }
                }
            }
        }
    }
    for (int thr = 10; thr<200; thr+=10){
        std::cout<<"Threshold "<<thr<<" mV:"<<std::endl;
        for (int ch = 0; ch<4; ch++){
            std::cout<<"  ch"<<'A'+ch<<" "<<counter[ch][thr]<<" "<<1.*counter[ch][thr]/nEntries/9.6<<" MHz"<<" occupancy "<<totalWidth[ch][thr]/nEntries/7600<<std::endl;
        }
    }

    return 0;
}
