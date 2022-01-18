#include "SvtPulseFitHistos.h"
#include <math.h>
#include "TCanvas.h"
#include "TRandom.h"

SvtPulseFitHistos::SvtPulseFitHistos(const std::string& inputName, ModuleMapper* mmapper) {
    m_name = inputName;
    mmapper_ = mmapper;
    svtid_map_ = mmapper_->buildChannelSvtIDMap();

    TH1F* hpe_h = new TH1F("nhits_per_event","nhits_per_event",160000,0,160000);
    histos1d_["nhits_per_event"] = hpe_h;

    TH2F* hpe_hh = new TH2F("nhits_per_event_hh","nhits_per_event_hh",160000,0,160000,8,-0.5,7.5);
    histos2d_["nhits_per_event_hh"] = hpe_hh;
    
}

SvtPulseFitHistos::~SvtPulseFitHistos() {
}

double SvtPulseFitHistos::GetHitTime(int sample_number, int cdel) {
    double time = 3.125*(8-cdel) + 24.0*sample_number;
    return time;
}

void SvtPulseFitHistos::buildRawSvtHitsTuple(std::vector<RawSvtHit*> *rawSvtHits_, FlatTupleMaker* rawhits_tup_) {

    /*
    bool debug = false;
    int nhits = rawSvtHits_->size();
    for (int i = 0; i < nhits; i++)
    {
        RawSvtHit* rawSvtHit = rawSvtHits_->at(i);
        int module = (rawSvtHit->getModule());
        int layer = (rawSvtHit->getLayer());
        std::string hwTag= mmapper_->getHwFromSw("ly"+std::to_string(layer)+"_m"+std::to_string(module));
        float channel = rawSvtHit->getStrip();
        int svtid = svtid_map_[hwTag].at(channel); 
        //int svtid = mmapper_->getSvtIDFromHWChannel(channel, hwTag, svtid_map_); 
        //int svtid = 99999;
        int calgroup = (int)channel%8;
        //Hard coded for now. Need to find way to determine this from data
        int cdel = 1;

        if(debug){
            std::cout << "module: " << module << " layer: " << layer << " channel: " << channel << 
                " svtid: " << svtid << " calgroup: " << calgroup << std::endl;
        }

        //ignore negative pulses
        if(rawSvtHit->getADCs()[3] < rawSvtHit->getADCs()[0])
            continue;

        rawhits_tup_->setVariableValue("hwTag", hwTag);
        rawhits_tup_->setVariableValue("layer", layer);
        rawhits_tup_->setVariableValue("module", module);
        rawhits_tup_->setVariableValue("channel", channel);
        rawhits_tup_->setVariableValue("svtid", svtid);
        rawhits_tup_->setVariableValue("calgroup", calgroup);
        rawhits_tup_->setVariableValue("cdel", cdel);

        rawhits_tup_->setVariableValue("adc0", rawSvtHit->getADCs()[0]);
        rawhits_tup_->setVariableValue("adc1", rawSvtHit->getADCs()[1]);
        rawhits_tup_->setVariableValue("adc2", rawSvtHit->getADCs()[2]);
        rawhits_tup_->setVariableValue("adc3", rawSvtHit->getADCs()[3]);
        rawhits_tup_->setVariableValue("adc4", rawSvtHit->getADCs()[4]);
        rawhits_tup_->setVariableValue("adc5", rawSvtHit->getADCs()[5]);

        rawhits_tup_->fill();
    }
    */
}

void SvtPulseFitHistos::defineTProfile(std::string name) {
    //TProfile* prof = new TProfile(name.c_str(),name.c_str(), 73,-1,145,2000,12000);
    //TProfile* prof = new TProfile(name.c_str(),name.c_str(), 48,-0.5,149,2000,12000);
    TProfile* prof = new TProfile(name.c_str(),name.c_str(), 49,-1.5625,151.5625,0,10000);
    tprofiles_[name] = prof;
}

void SvtPulseFitHistos::definePulseHistos(std::string name) {
    TH2F* prof = new TH2F(name.c_str(),name.c_str(), 49,-1.5625,151.5625,10000,0,10000);
    pulsehistos2d_[name] = prof;
}

void SvtPulseFitHistos::buildProfiles2019(TTree* rawhittree){

    std::cout << "WORKING ON 2019" << std::endl;
    double event, module, layer, channel, svtid, cdel, calgroup;
    double adc0, adc1, adc2, adc3, adc4, adc5;
    rawhittree->SetBranchAddress("event", &event);
    rawhittree->SetBranchAddress("svtid", &svtid);
    rawhittree->SetBranchAddress("channel", &channel);
    rawhittree->SetBranchAddress("layer", &layer);
    rawhittree->SetBranchAddress("module", &module);
    rawhittree->SetBranchAddress("adc0", &adc0);
    rawhittree->SetBranchAddress("adc1", &adc1);
    rawhittree->SetBranchAddress("adc2", &adc2);
    rawhittree->SetBranchAddress("adc3", &adc3);
    rawhittree->SetBranchAddress("adc4", &adc4);
    rawhittree->SetBranchAddress("adc5", &adc5);
    rawhittree->SetBranchAddress("cdel", &cdel);


    //std::map<int,std::pair<int,int>> cselmap = { {9,std::make_pair(0,2000)}, {8,std::make_pair(2000,4000)}, {7,std::make_pair(4000,6000)},{6,std::make_pair(6000,8000)},{5,std::make_pair(8000,10000)},{4,std::make_pair(10000,12000)},{3,std::make_pair(12000,14000)},{2,std::make_pair(14000,16000)},{1,std::make_pair(16000,18000)},{0,std::make_pair(18000,20000)} };
    //std::map<int,std::pair<int,int>> cselmap = {{0,std::make_pair(2000,4000)}, {1,std::make_pair(4000,6000)},{2,std::make_pair(6000,8000)},{3,std::make_pair(8000,10000)},{4,std::make_pair(10000,12000)},{5,std::make_pair(12000,14000)},{6,std::make_pair(14000,16000)},{7,std::make_pair(16000,18000)},{8,std::make_pair(18000,20000)} };
    //std::map<int,std::pair<int,int>> cselmap = {{0,std::make_pair(2000,4000)}, {1,std::make_pair(4000,6000)},{2,std::make_pair(6000,8000)},{3,std::make_pair(8000,10000)},{4,std::make_pair(10000,12000)},{5,std::make_pair(12000,14000)},{6,std::make_pair(14000,16000)},{7,std::make_pair(16000,18000)} };
    std::map<int,std::pair<int,int>> cgMap = { {0,std::make_pair(0,20000)},{1,std::make_pair(20000,40000)},{2,std::make_pair(40000,60000)},{3,std::make_pair(60000,80000)},{4,std::make_pair(80000,100000)},{5,std::make_pair(100000,120000)},{6,std::make_pair(120000,140000)},{7,std::make_pair(140000,160000)}  };
    std::map<int,std::pair<int,int>> cselmap = {{0,std::make_pair(0,2000)}, {1,std::make_pair(2000,4000)},{2,std::make_pair(4000,6000)},{3,std::make_pair(6000,8000)},{4,std::make_pair(8000,10000)},{5,std::make_pair(10000,12000)},{6,std::make_pair(12000,14000)},{7,std::make_pair(14000,16000)}, {8,std::make_pair(16000,18000)}, {-1, std::make_pair(18000,20000)} };

    int eventcount = -1;
    int onEvent = -1;
    int nrawhits = 0;
    int currentEvent = 0;

    long nentries = rawhittree->GetEntries();
    std::cout << "NENTRIES = " << nentries << std::endl;

    for(long i=0; i < nentries; i++){
        rawhittree->GetEntry(i);

        histos1d_["nhits_per_event"]->Fill(event);
        histos2d_["nhits_per_event_hh"]->Fill(event,(int)channel%8);
        continue;

        /*
        if(event > 0 && event < 10000){
            if((int)channel%8 == 1 && adc1 > 6500){
                std::cout << "event " << event << "| channel " << channel << " --> ADC0: " << adc0 << " ADC1: " << adc1 << " ADC2: " << adc2 << " ADC3: " << adc3 << " ADC4: " << adc4 << " ADC5: " << adc5 <<std::endl;
            }
        }
        else
            continue;
        */

        if(eventcount >= 20000){
            eventcount = -1;
        }

        if(event != onEvent){
            eventcount = eventcount + 2;
            onEvent = event;
        }

        //std::cout << "EVENTCOUNT: " << eventcount << std::endl;

        int modgroup = (int)channel%8;
        //std::cout << "modgroup: " << modgroup << std::endl;
        std::pair<int,int> cgRange = cgMap[modgroup];
        if(event >= cgRange.first && event < cgRange.second){
            bool inrange = true;
            //std::cout << "IS IN RANGE: " << cgRange.first << "-" << cgRange.second << std::endl;
        }
        else
            continue;

        std::map<int,std::pair<int,int>>::iterator it;
        int csel;
        for (it=cselmap.begin(); it != cselmap.end(); it++) {
            int min = it->second.first;
            int max = it->second.second;
            if (eventcount >= min && eventcount < max){
                csel = it->first;
                break;
            }
        }

        if(csel == -1 || csel == 0)
            continue;
        
        //std::cout << "EVENT " << event <<  "; event count: " << eventcount << std::endl;
        //std::cout << "CHANNEL : " << channel << std::endl;
        //std::cout << "channel%8: " << (int)channel%8 << std::endl;
        //if ((int)channel%8 != 0 && (int)channel%8 != 1)
        //    continue;

        //Build adcs vector
        std::vector<double> adcs;
        adcs.push_back(adc0);
        adcs.push_back(adc1);
        adcs.push_back(adc2);
        adcs.push_back(adc3);
        adcs.push_back(adc4);
        adcs.push_back(adc5);
        channel = (double)((int)channel);

        if (svtid > 1000)
            continue;
        //std::cout << "event : " << event << std::endl;

        std::string hwTag = mmapper_->getHwFromSw("ly"+std::to_string((int)layer)+"_m"+std::to_string((int)module));

        std::string name = hwTag+"_ch_"+std::to_string((int)channel)+"_svtid_"+std::to_string((int)svtid);
        //std::cout << name << std::endl;

        /*
        //check if tprofile exists
        if (tprofiles_.find(name) != tprofiles_.end()) {
            bool exists = true;
        }
        else{
           defineTProfile(name); 
        }
        */

        //Temporary analysis on 2d histos using incorrect but useful 25ns clock binning
        if (pulsehistos2d_.find(name) != pulsehistos2d_.end()){
            bool exists = true;
        }
        else
            definePulseHistos(name);


        for(int s=0; s < 6; s++){
            double adc = adcs.at(s);
            //double time = GetHitTime(s, csel);
            //tprofiles_[name]->Fill(time,adc,1.0);
            double time = 3.125*(8-csel) + 25.0*s;
            pulsehistos2d_[name]->Fill(time,adc);
        }
    }
}

void SvtPulseFitHistos::fitRawHitPulses(TTree* rawhittree, FlatTupleMaker* rawhitfits_tup) {

    double module, layer, channel, svtid, cdel, calgroup;
    double adc0, adc1, adc2, adc3, adc4, adc5;

    rawhittree->SetBranchAddress("svtid", &svtid);
    rawhittree->SetBranchAddress("channel", &channel);
    rawhittree->SetBranchAddress("layer", &layer);
    rawhittree->SetBranchAddress("module", &module);
    rawhittree->SetBranchAddress("adc0", &adc0);
    rawhittree->SetBranchAddress("adc1", &adc1);
    rawhittree->SetBranchAddress("adc2", &adc2);
    rawhittree->SetBranchAddress("adc3", &adc3);
    rawhittree->SetBranchAddress("adc4", &adc4);
    rawhittree->SetBranchAddress("adc5", &adc5);
    rawhittree->SetBranchAddress("cdel", &cdel);

    chi2_h_ = new TH1F("bestchi2","best_chi2;chi2;entries",1000,0,10);
    ndf_h_ = new TH1F("bestndf","best_ndf;ndf;entries",1000,0,10);
    t0_amp_h = new TH2F("t0_v_amp","t0_v_amp;t0 (ns);Amp (adc)",100,0,100,4000,0,400000);
    tau1_2_h = new TH2F("tau1_v_tau2","tau1_v_tau2;tau1;tau2",1000,0,100,1000,0,100);
    
    tau1_v_id = new TH1F("tau1_v_id","tau1_v_id",25000,-0.5,24999.5);
    tau2_v_id = new TH1F("tau2_v_id","tau2_v_id",25000,-0.5,24999.5);

    long nentries = rawhittree->GetEntries();
    for(long i=0; i < nentries; i++){
    //for(long i=0; i < 1000; i++){
        rawhittree->GetEntry(i);
        //Build adcs vector
        std::vector<double> adcs;
        adcs.push_back(adc0);
        adcs.push_back(adc1);
        adcs.push_back(adc2);
        adcs.push_back(adc3);
        adcs.push_back(adc4);
        adcs.push_back(adc5);
        channel = (double)((int)channel);

        std::string hwTag = mmapper_->getHwFromSw("ly"+std::to_string((int)layer)+"_m"+std::to_string((int)module));

        std::string name = hwTag+"_ch_"+std::to_string((int)channel)+"_svtid_"+std::to_string((int)svtid);
        //std::cout << name << std::endl;

        //check if tprofile exists
        if (tprofiles_.find(name) != tprofiles_.end()) {
            bool exists = true;
        }
        else{
           defineTProfile(name); 
        }

        //For 2021 Jlab cal pulse run, default cdel=0
        for(int t=0; t < 6; t++){
            double adc = adcs.at(t);
            double time = GetHitTime(t, 8);
            tprofiles_[name]->Fill(time,adc,1.0);
        }
    }

    for (ittp it = tprofiles_.begin(); it!=tprofiles_.end(); ++it) {
        /*
        std::string s = it->first;
        int svtid = std::stoi(s.substr(s.find("svtid_")+6,s.size()-1));
        std::string hwTag = s.substr(0,s.find("_"));
        int channel = std::stoi(s.substr(s.find("ch_"),s.find("_s")-5));
        */

        fitPulse(it->second, rawhitfits_tup);
    }
}

double SvtPulseFitHistos::getAmplitudeIntegralNorm(double tau1, double tau2) {
    double peak_t = 3.0*(pow((tau1* (pow(tau2,3))),(0.25)));
    double A = (pow(tau1,2))/(pow((tau1-tau2),3));
    double B = (tau1-tau2)/(tau1*tau2);
    double peakAmp = A * (exp(-peak_t / tau1) - exp(-peak_t/tau2) * (1 + peak_t * B + 0.5 * peak_t * peak_t *B*B));
    return peakAmp;
}

TF1* SvtPulseFitHistos::fourPoleFitFunction(){
    TF1* func = new TF1("pulsefit","(TMath::Max(x-[0],0.0)/(x-[0]))*([3])*(([1]^2)/(([1]-[2])^(3))) * ( exp(-(x-[0])/[1]) - ( exp(-(x-[0])/[2]) * ( ((((([1]-[2])/([1]*[2]))*(x-[0]))^(0))) +  ((((([1]-[2])/([1]*[2]))*(x-[0]))^(1))) + ((((([1]-[2])/([1]*[2]))*(x-[0]))^(2))/2) ) ) ) + [4]",0.0,150.0,"");
    return func;
}

void SvtPulseFitHistos::fitPulse(TProfile* tprofile, FlatTupleMaker* rawhitfits_tup_){

    std::string s = tprofile->GetName();
    int svtid = std::stoi(s.substr(s.find("svtid_")+6,s.size()-1));
    std::string hwTag = s.substr(0,s.find("_"));
    int channel = std::stoi(s.substr(s.find("ch_")+3,s.find("_s")-5));
    std::string sw = mmapper_->getSwFromHw(hwTag);
    int layer = std::stoi(sw.substr(sw.find("ly")+2,sw.find("_")-1));
    int module = std::stoi(sw.substr(sw.find("m")+1,sw.size()-1));

    double t0 = 0.0;
    double tau1 = 55.0;
    double tau2 = 8.0;
    double amp = 0.0;
    double baseline = 0.0;
    bool noPulse = false;

    //Get amplitude and baseline seed
    double maxamp = 0.0;
    double count = 0.0;
    for (int i = 0; i < tprofile->GetNbinsX(); i++) {
        double a = tprofile->GetBinContent(i+1);
        if(a > maxamp){
            maxamp = a;
        }

        if(i < 9){
            if(a > 0.0){
                baseline = baseline + a;
                count = count + 1.0;
            }
        }
    }

    //baseline = baseline/count;
    //Fix baseline to value of first bin
    baseline = tprofile->GetBinContent(1);
    amp = maxamp;

    //If tsample3 < tsample 2, no pulse in channel
    int fbin = tprofile->FindFirstBinAbove(0.0);
    /*
    std::cout << tprofile->GetBinContent(fbin) << std::endl;
    std::cout << tprofile->GetBinContent(fbin+12) << std::endl;
    std::cout << tprofile->GetBinContent(fbin+(24)) << std::endl;
    std::cout << tprofile->GetBinContent(fbin+(36)) << std::endl;
    std::cout << tprofile->GetBinContent(fbin+(48)) << std::endl;
    std::cout << tprofile->GetBinContent(fbin+(60)) << std::endl;
    */
    if ( (tprofile->GetBinContent(fbin+(2*6*3)) < tprofile->GetBinContent(fbin+(2*6*2)) ) && (tprofile->GetBinContent(fbin+(2*6*3)) < tprofile->GetBinContent(fbin+(2*6*4)) )) {
        noPulse = true;
        std::cout << "No Pulse in SVTID " << svtid << std::endl;
    }

    TF1* fitfunc = fourPoleFitFunction();

    TRandom *r1=new TRandom();

    //Fit with random seeds
    int not_nan_count = 0;
    bool isnan = false;
    double bestchi2 = 9999999.0;
    double besttau1 = tau1;
    double besttau2 = tau2;
    double bestt0 = t0;
    double bestamp = amp;
    double bestndf = 0.0;
    int maxiter = 20;
    int iter = 0;
    //while ((iter < maxiter || (not_nan_count < 5 && iter < 46)) || (isnan && iter < 40)){ 
    while(iter < maxiter){
        double rtau1;
        double rtau2;
        double rt0;
        double ramp; 

        if (svtid < 4096){
            rtau1 = r1->Gaus(53,2);
            rtau2 = r1->Gaus(8,2);
            rt0 = r1->Gaus(19,2);
            ramp = r1->Gaus(320000,2);
        }
        else {
            rtau1 = r1->Gaus(55,2);
            rtau2 = r1->Gaus(12,2);
            rt0 = r1->Gaus(16.5,2);
            ramp = r1->Gaus(300000,2);
        }

        fitfunc->SetParameter(1,rtau1);
        fitfunc->SetParameter(2,rtau2);
        fitfunc->SetParameter(0,rt0);
        fitfunc->SetParameter(3,ramp);
        fitfunc->FixParameter(4,baseline);

        tprofile->Fit(fitfunc, "q");

        //If fit gets nan errors, indicative of failed fit. 
        if (TMath::IsNaN(fitfunc->GetParError(0)) || TMath::IsNaN(fitfunc->GetParError(1)) || TMath::IsNaN(fitfunc->GetParError(2)) || TMath::IsNaN(fitfunc->GetParError(3))){
            isnan = true;
        }
        else{
            isnan = false;
            not_nan_count++;
        }

        double chi2 = fitfunc->GetChisquare();
        double ndf = (double) fitfunc->GetNDF();
        double chi2ndf = chi2/ndf;

        if(chi2ndf < bestchi2){
            bestchi2 = chi2ndf;
            besttau1 = rtau1;
            besttau2 = rtau2;
            bestt0 = rt0;
            bestamp = ramp;
            bestndf = ndf;
        }

        iter++;
        if(bestchi2 < 1)
            break;
    }

    fitfunc->SetParameter(1,besttau1);
    fitfunc->SetParameter(2,besttau2);
    fitfunc->SetParameter(0,bestt0);
    fitfunc->SetParameter(3,bestamp);
    fitfunc->FixParameter(4,baseline);
    tprofile->Fit(fitfunc, "q");

    t0 = fitfunc->GetParameter(0);
    tau1 = fitfunc->GetParameter(1);
    tau2 = fitfunc->GetParameter(2);
    amp = fitfunc->GetParameter(3);
    baseline = fitfunc->GetParameter(4);
    double chi2 = fitfunc->GetChisquare();
    double ndf = fitfunc->GetNDF();

    double t0err = fitfunc->GetParError(0);
    double tau1err = fitfunc->GetParError(1);
    double tau2err = fitfunc->GetParError(2);
    double amperr = fitfunc->GetParError(3);

    if (isnan){
        std::cout << "NAN FIT SVTID: " << svtid << std::endl;
    }
    if (TMath::IsNaN(fitfunc->GetParError(0)) || TMath::IsNaN(fitfunc->GetParError(1)) || TMath::IsNaN(fitfunc->GetParError(2)) || TMath::IsNaN(fitfunc->GetParError(3))){
        nan_channels_++;
    }

    rawhitfits_tup_->setVariableValue("hwTag", hwTag);
    rawhitfits_tup_->setVariableValue("channel", channel);
    rawhitfits_tup_->setVariableValue("svtid", svtid);
    rawhitfits_tup_->setVariableValue("module", module);
    rawhitfits_tup_->setVariableValue("layer", layer);
    rawhitfits_tup_->setVariableValue("noPulse", noPulse);
    rawhitfits_tup_->setVariableValue("t0", t0);
    rawhitfits_tup_->setVariableValue("tau1", tau1);
    rawhitfits_tup_->setVariableValue("tau2", tau2);
    rawhitfits_tup_->setVariableValue("amp", amp);
    rawhitfits_tup_->setVariableValue("baseline", baseline);
    rawhitfits_tup_->setVariableValue("chi2", chi2);
    rawhitfits_tup_->setVariableValue("ndf", ndf);
    rawhitfits_tup_->setVariableValue("t0err", t0err);
    rawhitfits_tup_->setVariableValue("tau1err", tau1err);
    rawhitfits_tup_->setVariableValue("tau2err", tau2err);
    rawhitfits_tup_->setVariableValue("amperr", amperr);
    rawhitfits_tup_->setVariableValue("integralNorm", SvtPulseFitHistos::getAmplitudeIntegralNorm(tau1, tau2));

    rawhitfits_tup_->fill();

    //Fill Histograms
    chi2_h_->Fill(fitfunc->GetChisquare()/fitfunc->GetNDF());
    t0_amp_h->Fill(fitfunc->GetParameter(0),fitfunc->GetParameter(3),1.0);
    tau1_2_h->Fill(fitfunc->GetParameter(1),fitfunc->GetParameter(2),1.0);
    tau1_v_id->SetBinContent(svtid+1,fitfunc->GetParameter(1));
    tau2_v_id->SetBinContent(svtid+1,fitfunc->GetParameter(2));

    tprofile->Draw();
}

void SvtPulseFitHistos::saveTProfiles(TFile* outF, std::string folder) {
    for (ittp it = tprofiles_.begin(); it!=tprofiles_.end(); ++it) {

        outF->cd();
        std::string s = it->first;
        std::string delim = "_ch";
        std::string dirname = s.substr(0,s.find(delim));

        bool exists = false;
        TKey* key;
        TIter next( outF->GetListOfKeys());
        while( (key = (TKey*) next())){
            std::string classname = key->GetClassName();
            if(classname.find("TDirectoryFile") != std::string::npos){
                std::string keyname = key->GetName();
                if(keyname.find(dirname) != std::string::npos){
                    exists = true;
                }
            }
        }
        //std::cout << outF->GetDirectory(dirname.c_str()) << std::endl;
        if(!exists){
            gDirectory->mkdir(dirname.c_str());
        }

        outF->cd(dirname.c_str());
        it->second->Write();
    }
}

void SvtPulseFitHistos::saveHistos(TFile* outFile){
    outFile->cd();

    std::cout << "Saving 2d histos" << std::endl;
    std::cout << "SIZE OF 2DHISTOS : " << pulsehistos2d_.size() << std::endl;
    for (itth it = pulsehistos2d_.begin(); it!=pulsehistos2d_.end(); ++it) {

        outFile->cd();
        std::string s = it->first;
        std::string delim = "_ch";
        std::string dirname = s.substr(0,s.find(delim));

        bool exists = false;
        TKey* key;
        TIter next( outFile->GetListOfKeys());
        while( (key = (TKey*) next())){
            std::string classname = key->GetClassName();
            if(classname.find("TDirectoryFile") != std::string::npos){
                std::string keyname = key->GetName();
                if(keyname.find(dirname) != std::string::npos){
                    exists = true;
                }
            }
        }
        //std::cout << outF->GetDirectory(dirname.c_str()) << std::endl;
        if(!exists){
            gDirectory->mkdir(dirname.c_str());
        }

        outFile->cd(dirname.c_str());
        it->second->Write();
    }

    //save 1d histos
    outFile->cd();
    typedef std::map<std::string, TH1F*>::iterator it1d;
    for (it1d it = histos1d_.begin(); it!=histos1d_.end(); it++) {
        it->second->Write();
    }
    //save 2d histos
    outFile->cd();
    typedef std::map<std::string, TH2F*>::iterator it2d;
    for (it2d it = histos2d_.begin(); it!=histos2d_.end(); it++) {
        it->second->Write();
    }

    /*
    chi2_h_->Write();
    ndf_h_->Write();
    t0_amp_h->Write();
    tau1_2_h->Write();
    tau1_v_id->Write();
    tau2_v_id->Write();
    */

    std::cout << "Number of NAN channel fits: " << nan_channels_ << std::endl;
}

