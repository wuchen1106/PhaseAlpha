#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <cstdlib>
#include <stdlib.h> /* atoi, atof */

#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLine.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>

#define NCH 6
int NCHmax = NCH;

int verbose = 0;

int iEntryStart = 0;
TString runname = "test";
double coincidenceWindowWidth = 5; // build event with coincidence +-5 ns
double tmin_trig = 360; // ns, trigger time no earlier
double tmax_trig = 400; // ns, trigger time no later
double t0[NCH] = {0,1,2,8,0,0}; // time offset to cancel cable delays
double hmin_trig[NCH] = {30,30,3000,3000,30,30}; // mV, threshold for trigger
double hmin[NCH] = {30,30,30,30,30,30}; // mV, threshold for event
TString triggerMode = "T0-T1-T2"; // pT0+T0-T1-T2, T0+T1+T2, pT0+T0+T1+T2

Int_t run = 0;
Int_t event = 0;
std::map<int,bool> used[NCH];
std::vector<double> * start[NCH] = {0}; // ns
std::vector<int> * width[NCH] = {0}; // number of ADC samples over threshold (set in wf-analyzer)
std::vector<int> * fwhm[NCH] = {0}; // number of ADC samples over half-height
std::vector<double> * height[NCH] = {0}; // mV over pedestal
std::vector<double> * area[NCH] = {0}; // pC over pedestal
int    nEvents = 0;
std::vector<double> * evt_t[NCH] = {0};
std::vector<double> * evt_h[NCH] = {0};
std::vector<int>    * evt_w[NCH] = {0};
std::vector<int>    * evt_f[NCH] = {0};
std::vector<double> * evt_a[NCH] = {0};
std::vector<double> * evt_t0mean = 0;
std::vector<double> * evt_h0mean = 0;
std::vector<double> * evt_w0mean = 0;
std::vector<double> * evt_f0mean = 0;
std::vector<double> * evt_a0mean = 0;
std::vector<int>    * evt_n0 = 0;
std::vector<double> * evt_tpmean = 0;
std::vector<double> * evt_hpmean = 0;
std::vector<double> * evt_wpmean = 0;
std::vector<double> * evt_fpmean = 0;
std::vector<double> * evt_apmean = 0;
std::vector<int>    * evt_np = 0;
std::vector<double> * evt_tmean = 0;
std::vector<double> * evt_hmean = 0;
std::vector<double> * evt_wmean = 0;
std::vector<double> * evt_fmean = 0;
std::vector<double> * evt_amean = 0;
std::vector<int>    * evt_n = 0;
std::vector<double> * evt_tmin = 0;
std::vector<double> * evt_tmax = 0;

bool removeEvent(int index){
    if (index<0||index>=nEvents) return true;
    nEvents--;
    for (int ch = 0; ch<NCHmax; ch++){
        evt_t[ch]->erase (evt_t[ch]->begin()+index);
        evt_h[ch]->erase (evt_h[ch]->begin()+index);
        evt_w[ch]->erase (evt_w[ch]->begin()+index);
        evt_f[ch]->erase (evt_f[ch]->begin()+index);
        evt_a[ch]->erase (evt_a[ch]->begin()+index);
    }
    evt_t0mean->erase (evt_t0mean->begin()+index);
    evt_h0mean->erase (evt_h0mean->begin()+index);
    evt_w0mean->erase (evt_w0mean->begin()+index);
    evt_f0mean->erase (evt_f0mean->begin()+index);
    evt_a0mean->erase (evt_a0mean->begin()+index);
    evt_n0->erase (evt_n0->begin()+index);
    evt_tpmean->erase (evt_tpmean->begin()+index);
    evt_hpmean->erase (evt_hpmean->begin()+index);
    evt_wpmean->erase (evt_wpmean->begin()+index);
    evt_fpmean->erase (evt_fpmean->begin()+index);
    evt_apmean->erase (evt_apmean->begin()+index);
    evt_np->erase (evt_np->begin()+index);
    evt_tmean->erase (evt_tmean->begin()+index);
    evt_hmean->erase (evt_hmean->begin()+index);
    evt_wmean->erase (evt_wmean->begin()+index);
    evt_fmean->erase (evt_fmean->begin()+index);
    evt_amean->erase (evt_amean->begin()+index);
    evt_n->erase (evt_n->begin()+index);
    evt_tmin->erase (evt_tmin->begin()+index);
    evt_tmax->erase (evt_tmax->begin()+index);
    return false;
};

struct Hit{
    int ch;
    int index;
    double t;
    double height;
    int    width;
    int    fwhm;
    double area;
    bool operator < (const Hit &b)const {
        return (t<b.t);
    };
    bool operator > (const Hit &b)const {
        return (t>b.t);
    };
    bool operator == (const Hit &b)const {
        return (t==b.t);
    };
    bool operator < (const Hit &b){
        return (t<b.t);
    };
    bool operator > (const Hit &b){
        return (t>b.t);
    };
    bool operator == (const Hit &b){
        return (t==b.t);
    };
    Hit(){
        ch = -1;
        index = -1;
        t = 0;
        height = 0;
        width = 0;
        fwhm = 0;
        area = 0;
    };
    void Print(){
        std::cout<<"ch "<<ch<<" ind. "<<index<<" H "<<height<<" t "<<t<<std::endl;
    };
};

/// A Combo is a cluster of hits coincide on time.
/// It may contain multiple hits in one channel, denoted by nHits[NCH]
/// The operator++ and isEnd flag is designed to loop all combinations by picking at most one hit per channel
struct Combo{
    int nHits[NCH];
    int iterator[NCH];
    bool isEnd;
    Combo& operator ++(){
        bool bringUp = true;
        for (int ch = 0; ch<NCHmax; ch++){
            if (nHits[ch]==0||iterator[ch]==nHits[ch]-1){
                bringUp = true;
                iterator[ch] = 0;
            }
            else{
                iterator[ch]++;
                break;
            }
            if (bringUp&&ch+1==NCHmax){
                isEnd = true;
            }
        }
        return *this;
    };
    Combo operator++(int)
    {
        Combo old = *this; // copy old value
        operator++();  // prefix increment
        return old;    // return old value
    };

    Combo(){
        for (int ch = 0; ch<NCHmax; ch++){
            nHits[ch] = 0; 
            iterator[ch] = 0;
        }
        isEnd = false;
    };
    void Print(){
        for (int ch = 0; ch<NCHmax; ch++){
            std::cout<<"ch "<<ch<<": "<<nHits[ch]<<" hits, iterator = "<<iterator[ch]<<std::endl;
        }
        std::cout<<"Is end? "<<(isEnd?"yes":"no")<<std::endl;
    };
};

struct Event{
    Hit hits[NCH];
    double tmin;
    double tmax;
    double averageH0;
    double averageW0;
    double averageF0;
    double averageA0;
    double averageT0; 
    int nHits0;
    double averageHp;
    double averageWp;
    double averageFp;
    double averageAp;
    double averageTp; 
    int nHitsp;
    double averageH;
    double averageW;
    double averageF;
    double averageA;
    double averageT; 
    int nHits;
    bool isValid;
    void update(){
        tmin = 1e9;
        for (int ch = 0; ch<NCHmax; ch++){
            if (hits[ch].height){
                if (hits[ch].t<tmin) tmin = hits[ch].t;
            }
        }
        tmax = -1e9;
        //for (int ch = 0; ch<NCHmax; ch++){
        //    if (hits[ch].height){
        //        if (hits[ch].t>tmin+coincidenceWindowWidth){
        //            if (verbose>5) std::cout<<"Hit removed! ch "<<ch<<" t "<<hits[ch].t<<">"<<tmin<<"+"<<coincidenceWindowWidth<<std::endl;
        //            Hit aHit;
        //            hits[ch] = aHit;
        //        }
        //        else if (hits[ch].t>tmax) tmax = hits[ch].t;
        //    }
        //}
        averageT0 = 0;
        averageH0 = 0;
        averageW0 = 0;
        averageF0 = 0;
        averageA0 = 0;
        nHits0 = 0;
        averageTp = 0;
        averageHp = 0;
        averageWp = 0;
        averageFp = 0;
        averageAp = 0;
        nHitsp = 0;
        averageT = 0;
        averageH = 0;
        averageW = 0;
        averageF = 0;
        averageA = 0;
        nHits = 0;
        for (int ch = 0; ch<NCHmax; ch++){
            if (ch<2){
                averageH0+=hits[ch].height;
                averageW0+=hits[ch].width;
                averageF0+=hits[ch].fwhm;
                averageA0+=hits[ch].area;
                if (hits[ch].height){
                    averageT0+=hits[ch].t;
                    nHits0++;
                }
            }
            else if (ch>=4){
                averageHp+=hits[ch].height;
                averageWp+=hits[ch].width;
                averageFp+=hits[ch].fwhm;
                averageAp+=hits[ch].area;
                if (hits[ch].height){
                    averageTp+=hits[ch].t;
                    nHitsp++;
                }
            }
            averageH+=hits[ch].height;
            averageW+=hits[ch].width;
            averageF+=hits[ch].fwhm;
            averageA+=hits[ch].area;
            if (hits[ch].height){
                averageT+=hits[ch].t;
                nHits++;
            }
        }
        if (nHits0>0){
            averageT0/=nHits0;
            averageH0/=nHits0;
            averageW0/=nHits0;
            averageF0/=nHits0;
            averageA0/=nHits0;
        }
        if (nHitsp>0){
            averageTp/=nHitsp;
            averageHp/=nHitsp;
            averageWp/=nHitsp;
            averageFp/=nHitsp;
            averageAp/=nHitsp;
        }
        if (nHits>0){
            averageT/=nHits;
            averageH/=nHits;
            averageW/=nHits;
            averageF/=nHits;
            averageA/=nHits;
        }
    };
    bool operator <(Event b){ // < means better
        update();
        b.update();
        if (nHits>b.nHits) return true;
        if (nHits<b.nHits) return false;
        if ((tmax-tmin)<=(b.tmax-b.tmin)&&(averageH0+averageHp+hits[2].height+hits[3].height)>=(b.averageH+b.averageHp+b.hits[2].height+b.hits[3].height)) return true;
        return (averageT<b.averageT);
    };
    Event(){
        averageT0 = 0;
        averageH0 = 0;
        averageW0 = 0;
        averageF0 = 0;
        averageA0 = 0;
        nHits0 = 0;
        averageTp = 0;
        averageHp = 0;
        averageWp = 0;
        averageFp = 0;
        averageAp = 0;
        nHitsp = 0;
        averageT = 0;
        averageH = 0;
        averageW = 0;
        averageF = 0;
        averageA = 0;
        nHits = 0;
    };
    void Print(){
        std::cout<<"Event with n "<<nHits<<" n0 "<<nHits0<<" np "<<nHitsp<<std::endl;
        for (int ch = 0; ch<NCHmax; ch++){
            hits[ch].Print();
        }
        std::cout<<"n0 = "<<nHits0<<" t0 = "<<averageT0<<", h0 = "<<averageH0<<std::endl;
        std::cout<<"np = "<<nHitsp<<" tp = "<<averageTp<<", hp = "<<averageHp<<std::endl;
    };
};

void SaveEvent(Event & event){
    for (int ch = 0; ch<NCHmax; ch++){
        if (event.hits[ch].index<0){
            evt_t[ch]->push_back(0);
            evt_h[ch]->push_back(0);
            evt_w[ch]->push_back(0);
            evt_f[ch]->push_back(0);
            evt_a[ch]->push_back(0);
        }
        else{
            int index = event.hits[ch].index;
            evt_t[ch]->push_back(start[ch]->at(index));
            evt_h[ch]->push_back(height[ch]->at(index));
            evt_w[ch]->push_back(width[ch]->at(index));
            evt_f[ch]->push_back(fwhm[ch]->at(index));
            evt_a[ch]->push_back(area[ch]->at(index));
        }
    }
    evt_t0mean->push_back(event.averageT0);
    evt_h0mean->push_back(event.averageH0);
    evt_w0mean->push_back(event.averageW0);
    evt_f0mean->push_back(event.averageF0);
    evt_a0mean->push_back(event.averageA0);
    evt_n0->push_back(event.nHits0);
    evt_tpmean->push_back(event.averageTp);
    evt_hpmean->push_back(event.averageHp);
    evt_wpmean->push_back(event.averageWp);
    evt_fpmean->push_back(event.averageFp);
    evt_apmean->push_back(event.averageAp);
    evt_np->push_back(event.nHitsp);
    evt_tmean->push_back(event.averageT);
    evt_hmean->push_back(event.averageH);
    evt_wmean->push_back(event.averageW);
    evt_fmean->push_back(event.averageF);
    evt_amean->push_back(event.averageA);
    evt_n->push_back(event.nHits);
    evt_tmin->push_back(event.tmin);
    evt_tmax->push_back(event.tmax);
    nEvents++;
}

std::vector<Hit> hitList;

bool analyzeCluster(int i_start, int i_stop){ // pick up according to candidate sorting
    bool ambiguous = false;
    std::vector<int> hitListIndexByChannel[NCH];
    for (int i = i_start; i<=i_stop; i++){
        int ch = hitList[i].ch;
        hitListIndexByChannel[ch].push_back(i);
    }
    int nChs = 0;
    int nHitsPerChMax = 0;
    for (int ch = 0; ch<NCHmax; ch++){
        if (hitListIndexByChannel[ch].size()) nChs++;
        if (hitListIndexByChannel[ch].size()>nHitsPerChMax) nHitsPerChMax = hitListIndexByChannel[ch].size();
    }
    if (nChs>1&&nHitsPerChMax>1) ambiguous = true;
    
    int nHits = i_stop-i_start+1;
    if (verbose>5){
        std::cout<<nHits<<" hits to be selected"<<std::endl;
        for (int ch = 0; ch<NCHmax; ch++){
            std::cout<<"ch "<<ch<<": "<<hitListIndexByChannel[ch].size()<<" hits"<<std::endl;
        }
    }
    while (nHits>0){
        if (verbose>5) std::cout<<" nHits = "<<nHits<<std::endl;
        Event theEvent;
        Combo combo;
        for (int ch = 0; ch<NCHmax; ch++){
            combo.nHits[ch] = hitListIndexByChannel[ch].size();
        }
        for (;!combo.isEnd;combo++){
            Event aEvent;
            for (int ch = 0; ch<NCHmax; ch++){
                int iterator = combo.iterator[ch];
                if (iterator<0||iterator>=combo.nHits[ch]) continue; // skip this channel; Possibly there is no hits in this channel
                int hitListIndex = hitListIndexByChannel[ch][iterator];
                int index = hitList[hitListIndex].index;
                if (!used[ch][index]){
                    aEvent.hits[ch] = hitList[hitListIndex];
                }
            }
            aEvent.update();
            if (verbose>5){
                for (int ch = 0; ch<NCHmax; ch++){
                    if (aEvent.hits[ch].height!=0){
                        std::cout<<"  Picked hit: ";
                        aEvent.hits[ch].Print();
                    }
                }
            }
            if (verbose>5) std::cout<<"Comparing theEvent and aEvent"<<std::endl;
            if (verbose>5) theEvent.Print();
            if (verbose>5) aEvent.Print();
            if (verbose>5) std::cout<<"theEvent empty?"<<(!theEvent.nHits?"yes":"no")<<" theEvent is worse?"<<(aEvent<theEvent?"yes":"no")<<std::endl;
            if (!theEvent.nHits||aEvent<theEvent){
                if (verbose>5) std::cout<<"Use this event"<<std::endl;
                theEvent = aEvent;
            }
        }
        //theEvent.Print();
        if (theEvent.nHits==0){
            std::cerr<<"ERROR: empty event!"<<std::endl;
            return false;
        }
        if (verbose>5) std::cout<<"  nHits "<<nHits;
        for (int ch = 0; ch<NCHmax; ch++){
            int index = theEvent.hits[ch].index;
            if (index<0) continue;
            used[ch][index] = true;
            nHits--;
        }
        if (verbose>5) std::cout<<" -> "<<nHits<<std::endl;
        SaveEvent(theEvent);
    }
    if (ambiguous) std::cout<<run<<" "<<event<<" "<<hitList[i_start].t<<" "<<hitList[i_stop].t<<std::endl;

    return ambiguous;
}

void print_usage(char * prog_name){
    fprintf(stderr,"Usage %s [options] runname\n",prog_name);
    fprintf(stderr,"[options]\n");
};

int main(int argc, char** argv)
{
    bool isOldData = false;

    // Load options
    int    opt_result;
    std::stringstream stream;
    while((opt_result=getopt(argc,argv,"oC:O:W:T:H:M:s:r:v:h"))!=-1){
        stream.clear();
        if(optarg) stream.str(optarg);
        switch(opt_result){
            /* INPUTS */
            case 'o':
                isOldData = true; // Februaray data with only 4 channels for picoscope: no pre-T0 yet
                break;
            case 's':
                iEntryStart = atoi(optarg);
                std::cerr<<"Starting from entry "<<iEntryStart<<std::endl;
                break;
            case 'C':
                coincidenceWindowWidth = atof(optarg);
                std::cerr<<"Coincidence window set to +-"<<coincidenceWindowWidth<<" ns"<<std::endl;
                break;
            case 'M':
                triggerMode = optarg;
                std::cerr<<"Trigger mode set to \""<<triggerMode<<"\""<<std::endl;
                break;
            case 'r':
                runname = optarg;
                std::cerr<<"Run name \""<<runname<<"\""<<std::endl;
                break;
            case 'W':
                stream>>tmin_trig>>tmax_trig;
                std::cerr<<"Set trigger window to "<<tmin_trig<<" ~ "<<tmax_trig<<" ns"<<std::endl;
                break;
            case 'H':
                std::cerr<<"Set peak hight cut for event building at";
                for (int i = 0; i<NCH; i++){
                    stream>>hmin[i];
                    std::cout<<" "<<hmin[i];
                }
                std::cerr<<" mV"<<std::endl;
                break;
            case 'O':
                std::cerr<<"Set time offset to each channel";
                for (int i = 0; i<NCH; i++){
                    stream>>t0[i];
                    std::cout<<" "<<t0[i];
                }
                std::cout<<" ns (to be subtracted)"<<std::endl;
                break;
            case 'T':
                std::cerr<<"Set peak hight cut for trigger at";
                for (int i = 0; i<NCH; i++){
                    stream>>hmin_trig[i];
                    std::cout<<" "<<hmin_trig[i];
                }
                std::cerr<<" mV"<<std::endl;
                break;
            case 'v':
                verbose = atoi(optarg);
                std::cerr<<"verbose level set to "<<verbose<<std::endl;
                break;
            case '?':
                fprintf(stderr,"Wrong option! optopt=%c, optarg=%s\n", optopt, optarg);
            case 'h':
            default:
                print_usage(argv[0]);
                return 1;
        }
    }
    if (isOldData) NCHmax = 4;

    TChain* chain = new TChain("t");
    for (int i = optind; i<argc; i++){
        chain->Add(argv[i]);
    }
    if(0 == chain->GetNtrees()) return 1;

    // Set Branch Addresses
    double rms[NCH] = {0};
    double pedestal[NCH] = {0};
    int    nPeaks[NCH] = {0};
    chain->SetBranchAddress("run",&run);
    chain->SetBranchAddress("event",&event);
    for (int ch = 0; ch<NCHmax; ch++){
        chain->SetBranchAddress(Form("rms%c",'A'+ch),&rms[ch]);
        chain->SetBranchAddress(Form("ped%c",'A'+ch),&pedestal[ch]);
        chain->SetBranchAddress(Form("nPeaks%c",'A'+ch),&nPeaks[ch]);
        chain->SetBranchAddress(Form("start%c",'A'+ch),&start[ch]);
        chain->SetBranchAddress(Form("width%c",'A'+ch),&width[ch]);
        chain->SetBranchAddress(Form("fwhm%c",'A'+ch),&fwhm[ch]);
        chain->SetBranchAddress(Form("height%c",'A'+ch),&height[ch]);
        chain->SetBranchAddress(Form("area%c",'A'+ch),&area[ch]);
    }

    TFile * ofile = new TFile(Form("anaout/%s.event.root",runname.Data()),"recreate");
    TTree * otree = new TTree("t","t");
    otree->Branch("run",&run);
    otree->Branch("event",&event);
    int foundTrigger = 0;
    double trig_t[NCH] = {0};
    double trig_h[NCH] = {0};
    int    trig_w[NCH] = {0};
    int    trig_f[NCH] = {0};
    double trig_a[NCH] = {0};
    double trig_t0mean = 0;
    double trig_h0mean = 0;
    double trig_w0mean = 0;
    double trig_f0mean = 0;
    double trig_a0mean = 0;
    int    trig_n0     = 0;
    double trig_tpmean = 0;
    double trig_hpmean = 0;
    double trig_wpmean = 0;
    double trig_fpmean = 0;
    double trig_apmean = 0;
    int    trig_np     = 0;
    double trig_tmean = 0;
    double trig_hmean = 0;
    double trig_wmean = 0;
    double trig_fmean = 0;
    double trig_amean = 0;
    int    trig_n     = 0;
    bool foundPeak110 = false; // peak at 111 ns
    bool foundPeak470 = false; // peak at 470 ns
    double peak110_t[NCH] = {0};
    double peak110_h[NCH] = {0};
    int    peak110_w[NCH] = {0};
    int    peak110_f[NCH] = {0};
    double peak110_a[NCH] = {0};
    double peak110_t0mean = 0;
    double peak110_h0mean = 0;
    double peak110_w0mean = 0;
    double peak110_f0mean = 0;
    double peak110_a0mean = 0;
    int    peak110_n0    = 0;
    double peak110_tpmean = 0;
    double peak110_hpmean = 0;
    double peak110_wpmean = 0;
    double peak110_fpmean = 0;
    double peak110_apmean = 0;
    int    peak110_np    = 0;
    double peak110_tmean = 0;
    double peak110_hmean = 0;
    double peak110_wmean = 0;
    double peak110_fmean = 0;
    double peak110_amean = 0;
    int    peak110_n    = 0;
    double peak470_t[NCH] = {0};
    double peak470_h[NCH] = {0};
    int    peak470_w[NCH] = {0};
    int    peak470_f[NCH] = {0};
    double peak470_a[NCH] = {0};
    double peak470_t0mean = 0;
    double peak470_h0mean = 0;
    double peak470_w0mean = 0;
    double peak470_f0mean = 0;
    double peak470_a0mean = 0;
    int    peak470_n0    = 0;
    double peak470_tpmean = 0;
    double peak470_hpmean = 0;
    double peak470_wpmean = 0;
    double peak470_fpmean = 0;
    double peak470_apmean = 0;
    int    peak470_np    = 0;
    double peak470_tmean = 0;
    double peak470_hmean = 0;
    double peak470_wmean = 0;
    double peak470_fmean = 0;
    double peak470_amean = 0;
    int    peak470_n    = 0;
    otree->Branch("nEvents",&nEvents);
    for (int ch = 0; ch<NCHmax; ch++){
        otree->Branch(Form("trig_t%c",'A'+ch),&trig_t[ch]);
        otree->Branch(Form("trig_h%c",'A'+ch),&trig_h[ch]);
        otree->Branch(Form("trig_w%c",'A'+ch),&trig_w[ch]);
        otree->Branch(Form("trig_f%c",'A'+ch),&trig_f[ch]);
        otree->Branch(Form("trig_a%c",'A'+ch),&trig_a[ch]);
        otree->Branch(Form("peak110_t%c",'A'+ch),&peak110_t[ch]);
        otree->Branch(Form("peak110_h%c",'A'+ch),&peak110_h[ch]);
        otree->Branch(Form("peak110_w%c",'A'+ch),&peak110_w[ch]);
        otree->Branch(Form("peak110_f%c",'A'+ch),&peak110_f[ch]);
        otree->Branch(Form("peak110_a%c",'A'+ch),&peak110_a[ch]);
        otree->Branch(Form("peak470_t%c",'A'+ch),&peak470_t[ch]);
        otree->Branch(Form("peak470_h%c",'A'+ch),&peak470_h[ch]);
        otree->Branch(Form("peak470_w%c",'A'+ch),&peak470_w[ch]);
        otree->Branch(Form("peak470_f%c",'A'+ch),&peak470_f[ch]);
        otree->Branch(Form("peak470_a%c",'A'+ch),&peak470_a[ch]);
        otree->Branch(Form("evt_t%c",'A'+ch),&evt_t[ch]);
        otree->Branch(Form("evt_h%c",'A'+ch),&evt_h[ch]);
        otree->Branch(Form("evt_w%c",'A'+ch),&evt_w[ch]);
        otree->Branch(Form("evt_f%c",'A'+ch),&evt_f[ch]);
        otree->Branch(Form("evt_a%c",'A'+ch),&evt_a[ch]);
        // can skip the following hit level informations
        otree->Branch(Form("rms%c",'A'+ch),&rms[ch]);
        otree->Branch(Form("ped%c",'A'+ch),&pedestal[ch]);
        evt_t[ch] = new std::vector<double>;
        evt_h[ch] = new std::vector<double>;
        evt_w[ch] = new std::vector<int>;
        evt_f[ch] = new std::vector<int>;
        evt_a[ch] = new std::vector<double>;
    }
    otree->Branch("evt_t0",&evt_t0mean);
    otree->Branch("evt_h0",&evt_h0mean);
    otree->Branch("evt_w0",&evt_w0mean);
    otree->Branch("evt_f0",&evt_f0mean);
    otree->Branch("evt_a0",&evt_a0mean);
    otree->Branch("evt_n0",&evt_n0);
    otree->Branch("evt_tp",&evt_tpmean);
    otree->Branch("evt_hp",&evt_hpmean);
    otree->Branch("evt_wp",&evt_wpmean);
    otree->Branch("evt_fp",&evt_fpmean);
    otree->Branch("evt_ap",&evt_apmean);
    otree->Branch("evt_np",&evt_np);
    otree->Branch("evt_t",&evt_tmean);
    otree->Branch("evt_h",&evt_hmean);
    otree->Branch("evt_w",&evt_wmean);
    otree->Branch("evt_f",&evt_fmean);
    otree->Branch("evt_a",&evt_amean);
    otree->Branch("evt_n",&evt_n);
    otree->Branch("evt_tmin",&evt_tmin);
    otree->Branch("evt_tmax",&evt_tmax);
    evt_t0mean = new std::vector<double>;
    evt_h0mean = new std::vector<double>;
    evt_w0mean = new std::vector<double>;
    evt_f0mean = new std::vector<double>;
    evt_a0mean = new std::vector<double>;
    evt_n0 = new std::vector<int>;
    evt_tpmean = new std::vector<double>;
    evt_hpmean = new std::vector<double>;
    evt_wpmean = new std::vector<double>;
    evt_fpmean = new std::vector<double>;
    evt_apmean = new std::vector<double>;
    evt_np = new std::vector<int>;
    evt_tmean = new std::vector<double>;
    evt_hmean = new std::vector<double>;
    evt_wmean = new std::vector<double>;
    evt_fmean = new std::vector<double>;
    evt_amean = new std::vector<double>;
    evt_n = new std::vector<int>;
    evt_tmin = new std::vector<double>;
    evt_tmax = new std::vector<double>;
    otree->Branch("foundTrigger",&foundTrigger);
    otree->Branch("foundPeak110",&foundPeak110);
    otree->Branch("foundPeak470",&foundPeak470);
    otree->Branch("trig_t0",&trig_t0mean);
    otree->Branch("trig_h0",&trig_h0mean);
    otree->Branch("trig_w0",&trig_w0mean);
    otree->Branch("trig_f0",&trig_f0mean);
    otree->Branch("trig_a0",&trig_a0mean);
    otree->Branch("trig_n0",&trig_n0);
    otree->Branch("trig_tp",&trig_tpmean);
    otree->Branch("trig_hp",&trig_hpmean);
    otree->Branch("trig_wp",&trig_wpmean);
    otree->Branch("trig_fp",&trig_fpmean);
    otree->Branch("trig_ap",&trig_apmean);
    otree->Branch("trig_np",&trig_np);
    otree->Branch("trig_t",&trig_tmean);
    otree->Branch("trig_h",&trig_hmean);
    otree->Branch("trig_w",&trig_wmean);
    otree->Branch("trig_f",&trig_fmean);
    otree->Branch("trig_a",&trig_amean);
    otree->Branch("trig_n",&trig_n);
    otree->Branch("peak110_t0",&peak110_t0mean);
    otree->Branch("peak110_h0",&peak110_h0mean);
    otree->Branch("peak110_w0",&peak110_w0mean);
    otree->Branch("peak110_f0",&peak110_f0mean);
    otree->Branch("peak110_a0",&peak110_a0mean);
    otree->Branch("peak110_n0",&peak110_n0);
    otree->Branch("peak110_tp",&peak110_tpmean);
    otree->Branch("peak110_hp",&peak110_hpmean);
    otree->Branch("peak110_wp",&peak110_wpmean);
    otree->Branch("peak110_fp",&peak110_fpmean);
    otree->Branch("peak110_ap",&peak110_apmean);
    otree->Branch("peak110_np",&peak110_np);
    otree->Branch("peak110_t",&peak110_tmean);
    otree->Branch("peak110_h",&peak110_hmean);
    otree->Branch("peak110_w",&peak110_wmean);
    otree->Branch("peak110_f",&peak110_fmean);
    otree->Branch("peak110_a",&peak110_amean);
    otree->Branch("peak110_n",&peak110_n);
    otree->Branch("peak470_t0",&peak470_t0mean);
    otree->Branch("peak470_h0",&peak470_h0mean);
    otree->Branch("peak470_w0",&peak470_w0mean);
    otree->Branch("peak470_f0",&peak470_f0mean);
    otree->Branch("peak470_a0",&peak470_a0mean);
    otree->Branch("peak470_n0",&peak470_n0);
    otree->Branch("peak470_tp",&peak470_tpmean);
    otree->Branch("peak470_hp",&peak470_hpmean);
    otree->Branch("peak470_wp",&peak470_wpmean);
    otree->Branch("peak470_fp",&peak470_fpmean);
    otree->Branch("peak470_ap",&peak470_apmean);
    otree->Branch("peak470_np",&peak470_np);
    otree->Branch("peak470_t",&peak470_tmean);
    otree->Branch("peak470_h",&peak470_hmean);
    otree->Branch("peak470_w",&peak470_wmean);
    otree->Branch("peak470_f",&peak470_fmean);
    otree->Branch("peak470_a",&peak470_amean);
    otree->Branch("peak470_n",&peak470_n);

    Long64_t nEntries = chain->GetEntries();
    int nEntriesWithoutTrigger = 0;
    int nEntriesWithMultipleTrigger = 0;
    for(Long64_t entry=iEntryStart; entry<nEntries; entry++){
        chain->GetEntry(entry);
        if (verbose>1) std::cout<<"############ entry "<<entry<<" ###############"<<std::endl;

        nEvents = 0;
        foundTrigger = 0;
        foundPeak110 = false;
        foundPeak470 = false;
        trig_t0mean = 0;
        trig_h0mean = 0;
        trig_w0mean = 0;
        trig_f0mean = 0;
        trig_n0 = 0;
        trig_tpmean = 0;
        trig_hpmean = 0;
        trig_wpmean = 0;
        trig_fpmean = 0;
        trig_np = 0;
        trig_tmean = 0;
        trig_hmean = 0;
        trig_wmean = 0;
        trig_fmean = 0;
        trig_n = 0;
        peak110_t0mean = 0;
        peak110_h0mean = 0;
        peak110_w0mean = 0;
        peak110_f0mean = 0;
        peak110_n0 = 0;
        peak110_tpmean = 0;
        peak110_hpmean = 0;
        peak110_wpmean = 0;
        peak110_fpmean = 0;
        peak110_np = 0;
        peak110_tmean = 0;
        peak110_hmean = 0;
        peak110_wmean = 0;
        peak110_fmean = 0;
        peak110_n = 0;
        peak470_t0mean = 0;
        peak470_h0mean = 0;
        peak470_w0mean = 0;
        peak470_f0mean = 0;
        peak470_n0 = 0;
        peak470_tpmean = 0;
        peak470_hpmean = 0;
        peak470_wpmean = 0;
        peak470_fpmean = 0;
        peak470_np = 0;
        peak470_tmean = 0;
        peak470_hmean = 0;
        peak470_wmean = 0;
        peak470_fmean = 0;
        peak470_n = 0;
        for (int ch = 0; ch<NCHmax; ch++){
            trig_t[ch] = 0;
            trig_h[ch] = 0;
            trig_w[ch] = 0;
            trig_f[ch] = 0;
            trig_a[ch] = 0;
            peak110_t[ch] = 0;
            peak110_h[ch] = 0;
            peak110_w[ch] = 0;
            peak110_f[ch] = 0;
            peak110_a[ch] = 0;
            peak470_t[ch] = 0;
            peak470_h[ch] = 0;
            peak470_w[ch] = 0;
            peak470_f[ch] = 0;
            peak470_a[ch] = 0;
            evt_t[ch]->clear();
            evt_h[ch]->clear();
            evt_w[ch]->clear();
            evt_f[ch]->clear();
            evt_a[ch]->clear();
            used[ch].clear();
            hitList.clear();
            // add time offset
            for (unsigned int i = 0; i<start[ch]->size(); i++){
                start[ch]->at(i) -= t0[ch];
            }
        }
        evt_t0mean->clear();
        evt_h0mean->clear();
        evt_w0mean->clear();
        evt_f0mean->clear();
        evt_a0mean->clear();
        evt_n0->clear();
        evt_tpmean->clear();
        evt_hpmean->clear();
        evt_wpmean->clear();
        evt_fpmean->clear();
        evt_apmean->clear();
        evt_np->clear();
        evt_tmean->clear();
        evt_hmean->clear();
        evt_wmean->clear();
        evt_fmean->clear();
        evt_amean->clear();
        evt_n->clear();
        evt_tmin->clear();
        evt_tmax->clear();

        // form the events
        for (int ch = 0; ch<NCHmax; ch++){
            for (unsigned int i = 0; i<start[ch]->size(); i++){
                if (height[ch]->at(i)<hmin[ch]) continue;
                Hit aHit;
                aHit.ch = ch;
                aHit.height = height[ch]->at(i);
                aHit.width  = width[ch]->at(i);
                aHit.fwhm   = fwhm[ch]->at(i);
                aHit.area   = area[ch]->at(i);
                aHit.index = i;
                aHit.t = start[ch]->at(i);
                hitList.push_back(aHit);
            }
        }

        std::sort(hitList.begin(),hitList.end());

        double t_prev = -1e9;
        int i_start = -1;
        for (unsigned int i = 0; i<hitList.size(); i++){
            if (hitList[i].t-t_prev>coincidenceWindowWidth){
                if (i_start!=-1){ // analyze the current cluster
                    int i_stop = i-1;
                    analyzeCluster(i_start,i_stop);
                }
                i_start = i;
            }
            t_prev = hitList[i].t;
        }
        if (i_start!=-1){
            int i_stop = hitList.size()-1;
            analyzeCluster(i_start,i_stop);
        }

        // check which is the trigger
        int itrig = -1;
        for (int i = 0; i<nEvents; i++){
            bool isTrigger = true;
            if (evt_n0->at(i)<=0) isTrigger = false;
            if (evt_t0mean->at(i)<tmin_trig||evt_t0mean->at(i)>tmax_trig){
                isTrigger = false;
            }
            if (evt_h[0]->at(i)<hmin_trig[0]||evt_h[1]->at(i)<hmin_trig[1]){
                isTrigger = false;
            }
            if (!isOldData&&(triggerMode=="pT0+T0-T1-T2"||triggerMode=="pT0+T0+T1+T2")){
                if (evt_np->at(i)<=0) isTrigger = false;
                if (evt_tpmean->at(i)<tmin_trig||evt_tpmean->at(i)>tmax_trig){
                    isTrigger = false;
                }
                if (evt_h[4]->at(i)<hmin_trig[4]||evt_h[5]->at(i)<hmin_trig[5]){
                    isTrigger = false;
                }
            }
            if (triggerMode=="T0-T1-T2"||triggerMode=="pT0+T0-T1-T2"){ // require
                if (evt_h[2]->at(i)>hmin_trig[2]||evt_h[3]->at(i)>hmin_trig[3]){
                    isTrigger = false;
                }
            }
            else{
                if (evt_h[2]->at(i)<hmin_trig[2]||evt_h[3]->at(i)<hmin_trig[3]){
                    isTrigger = false;
                }
            }
            if (isTrigger){
                // update?
                if (foundTrigger==0){
                    bool update = false;
                    if (triggerMode=="T0-T1-T2"){
                        if (trig_n0<evt_n0->at(i)||trig_h0mean<evt_h0mean->at(i)){
                            update = true;
                        }
                    }
                    else if (triggerMode=="pT0+T0-T1-T2"){
                        if (trig_n0+trig_np<evt_n0->at(i)+evt_np->at(i)||trig_h0mean+trig_hpmean<evt_h0mean->at(i)+evt_hpmean->at(i)){
                            update = true;
                        }
                    }
                    else{
                        if (trig_n<evt_n->at(i)||trig_hmean<evt_hmean->at(i)){ // the new one is closer to tmax_trig
                            update = true;
                        }
                    }
                    if (update){
                        for (int ch = 0; ch<NCHmax; ch++){
                            trig_t[ch] = evt_t[ch]->at(i);
                            trig_h[ch] = evt_h[ch]->at(i);
                            trig_w[ch] = evt_w[ch]->at(i);
                            trig_f[ch] = evt_f[ch]->at(i);
                            trig_a[ch] = evt_a[ch]->at(i);
                        }
                        trig_t0mean = evt_t0mean->at(i);
                        trig_h0mean = evt_h0mean->at(i);
                        trig_w0mean = evt_w0mean->at(i);
                        trig_f0mean = evt_f0mean->at(i);
                        trig_a0mean = evt_a0mean->at(i);
                        trig_n0 = evt_n0->at(i);
                        trig_tpmean = evt_tpmean->at(i);
                        trig_hpmean = evt_hpmean->at(i);
                        trig_wpmean = evt_wpmean->at(i);
                        trig_fpmean = evt_fpmean->at(i);
                        trig_apmean = evt_apmean->at(i);
                        trig_np = evt_np->at(i);
                        trig_tmean = evt_tmean->at(i);
                        trig_hmean = evt_hmean->at(i);
                        trig_wmean = evt_wmean->at(i);
                        trig_fmean = evt_fmean->at(i);
                        trig_amean = evt_amean->at(i);
                        trig_n = evt_n->at(i);
                        itrig=i;
                    }
                }
                foundTrigger++;
            }
        }
        // check if we have the special peaks
        int i110 = -1;
        int i470 = -1;
        for (int i = 0; i<nEvents; i++){
            for (int ch = 0; ch<=1; ch++){
                if (evt_h[ch]->at(i)<=0) continue;
                double deltaT = evt_t[ch]->at(i)-trig_tmean;
                if (deltaT>=111&&deltaT<=116){
                    foundPeak110 = true;
                    i110 = i;
                }
                if (deltaT>=471&&deltaT<=476){
                    foundPeak470 = true;
                    i470 = i;
                }
            }
        }
        if (foundPeak110){
            for (int ch = 0; ch<NCHmax; ch++){
                peak110_t[ch] = evt_t[ch]->at(i110);
                peak110_h[ch] = evt_h[ch]->at(i110);
                peak110_w[ch] = evt_w[ch]->at(i110);
                peak110_f[ch] = evt_f[ch]->at(i110);
                peak110_a[ch] = evt_a[ch]->at(i110);
            }
            peak110_t0mean = evt_t0mean->at(i110);
            peak110_h0mean = evt_h0mean->at(i110);
            peak110_w0mean = evt_w0mean->at(i110);
            peak110_f0mean = evt_f0mean->at(i110);
            peak110_a0mean = evt_a0mean->at(i110);
            peak110_n0 = evt_n0->at(i110);
            peak110_tpmean = evt_tpmean->at(i110);
            peak110_hpmean = evt_hpmean->at(i110);
            peak110_wpmean = evt_wpmean->at(i110);
            peak110_fpmean = evt_fpmean->at(i110);
            peak110_apmean = evt_apmean->at(i110);
            peak110_np = evt_np->at(i110);
            peak110_tmean = evt_tmean->at(i110);
            peak110_hmean = evt_hmean->at(i110);
            peak110_wmean = evt_wmean->at(i110);
            peak110_fmean = evt_fmean->at(i110);
            peak110_amean = evt_amean->at(i110);
            peak110_n = evt_n->at(i110);
        }
        if (foundPeak470){
            for (int ch = 0; ch<NCHmax; ch++){
                peak470_t[ch] = evt_t[ch]->at(i470);
                peak470_h[ch] = evt_h[ch]->at(i470);
                peak470_w[ch] = evt_w[ch]->at(i470);
                peak470_f[ch] = evt_f[ch]->at(i470);
                peak470_a[ch] = evt_a[ch]->at(i470);
            }
            peak470_t0mean = evt_t0mean->at(i470);
            peak470_h0mean = evt_h0mean->at(i470);
            peak470_w0mean = evt_w0mean->at(i470);
            peak470_f0mean = evt_f0mean->at(i470);
            peak470_a0mean = evt_a0mean->at(i470);
            peak470_n0 = evt_n0->at(i470);
            peak470_tpmean = evt_tpmean->at(i470);
            peak470_hpmean = evt_hpmean->at(i470);
            peak470_wpmean = evt_wpmean->at(i470);
            peak470_fpmean = evt_fpmean->at(i470);
            peak470_apmean = evt_apmean->at(i470);
            peak470_np = evt_np->at(i470);
            peak470_tmean = evt_tmean->at(i470);
            peak470_hmean = evt_hmean->at(i470);
            peak470_wmean = evt_wmean->at(i470);
            peak470_fmean = evt_fmean->at(i470);
            peak470_amean = evt_amean->at(i470);
            peak470_n = evt_n->at(i470);
        }

        /*
        std::vector<int> indexToRemove;
        if(itrig>=0) indexToRemove.push_back(itrig);
        if(i110>=0) indexToRemove.push_back(i110);
        if(i470>=0) indexToRemove.push_back(i470);
        // now sort them from large to small, so that we can remove the indices one by one without breaking the rest indices
        std::sort(indexToRemove.begin(),indexToRemove.end(),std::greater<int>());
        for (unsigned int ir = 0; ir<indexToRemove.size(); ir++){
            if (removeEvent(indexToRemove[ir])){
                std::cerr<<"Error while removing event!"<<std::endl;
                return -1;
            }
        }
        */

        if (foundTrigger==0) nEntriesWithoutTrigger++;
        else if (foundTrigger>1) nEntriesWithMultipleTrigger++;

        otree->Fill();
    }

    std::cout<<"  "<<nEntries<<" entries, "<<nEntriesWithoutTrigger<<" with no trigger "<<nEntriesWithMultipleTrigger<<" with more than one trigger"<<std::endl;

    otree->Write();
    ofile->Close();

    return 0;
}
