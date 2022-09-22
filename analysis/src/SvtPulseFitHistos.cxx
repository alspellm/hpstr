#include "SvtPulseFitHistos.h"
#include <math.h>
#include "TCanvas.h"
#include "TRandom.h"

SvtPulseFitHistos::SvtPulseFitHistos(const std::string& inputName, ModuleMapper* mmapper, int year) {
    m_name = inputName;
    mmapper_ = mmapper;
    svtid_map_ = mmapper_->buildChannelSvtIDMap();
    year_ = year;
}

SvtPulseFitHistos::~SvtPulseFitHistos() {
}

void SvtPulseFitHistos::initHistos(){
    //1d histos
    TH1F* hpe_h = new TH1F("nhits_per_event","nhits_per_event",25000,0,25000);
    TH1F* chi2_h = new TH1F("bestchi2ndf","best_chi2;chi2;entries",2000,0,20);
    TH1F* tau1_v_id = new TH1F("tau1_v_id","tau1_v_id",25000,-0.5,24999.5);
    TH1F* tau2_v_id = new TH1F("tau2_v_id","tau2_v_id",25000,-0.5,24999.5);

    histos1d_["nhits_per_event"] = hpe_h;
    histos1d_["chi2"] = chi2_h;
    histos1d_["tau1_v_id"] = tau1_v_id;
    histos1d_["tau2_v_id"] = tau2_v_id;

    //2d histos
    TH2F* adcs_hh = new TH2F("adcs","adcs;sample N;ADC",6,0,6,5000,0,10000);
    TH2F* t0_v_amp_hh = new TH2F("t0_v_amp","t0_v_amp;t0 (ns);Amp (adc)",300,-20,10,4000,0,400000);
    TH2F* tau1_v_tau2_hh = new TH2F("tau1_v_tau2","tau1_v_tau2;tau1;tau2",1000,0,100,1000,0,100);

    histos2d_["adcs"] = adcs_hh;
    histos2d_["t0_v_amp"] = t0_v_amp_hh;
    histos2d_["tau1_v_tau2"] = tau1_v_tau2_hh;
}

double SvtPulseFitHistos::GetHitTime(int sample_number, int cdel) {
    double time = 3.125*(8-cdel) + 24.0*sample_number;
    return time;
}

void SvtPulseFitHistos::definePulseHistos(std::string name) {

    //Binned for the 24ns JLab Clock with 3 ns APV25 internal delay 
    TH2F* prof = new TH2F(name.c_str(),name.c_str(), 51,-1.5,151.5,10000,0,10000);
    pulsehistos2d_[name] = prof;


    //store channel baselines in histogram
    TH1F* baseline = new TH1F((name+"_baseline").c_str(),name.c_str(),10000,0,10000);
    baselines_[name] = baseline;
}

void SvtPulseFitHistos::buildPulsesFromTree(TTree* rawhitsTree){
    std::cout << "Building pulses for  year " << year_ << std::endl;
    if(year_ == 2019){
        jlab2019CalPulseScan(rawhitsTree);
    }
    else if(year_ == 2021){
        jlab2021CalPulseScan(rawhitsTree);
    }
    else {
        std::cout << "ERROR! YEAR " << year_ << "DOES NOT HAVE METHOD TO PARSE PULSE DATA AVAILABLE. CLOSING" << std::endl;
    }
}

void SvtPulseFitHistos::jlab2021CalPulseScan(TTree* rawhitsTree){
    double module, layer, channel, svtid, cdel, calgroup;
    double adc0, adc1, adc2, adc3, adc4, adc5;

    //Read rawsvthit TTree
    rawhitsTree->Print();
    std::vector<RawSvtHit*>* rawSvtHits{};
    int event;
    std::cout << "SETTING BRANCHES TO READ" << std::endl;
    rawhitsTree->SetBranchAddress("event", &event);
    rawhitsTree->SetBranchAddress("rawsvthits", &rawSvtHits);
    std::cout << "BRANCHES READ" << std::endl;

    int nentries = rawhitsTree->GetEntries();
    std::cout << "number events: " << nentries;
    for(int i=0; i < nentries; i++){
        rawhitsTree->GetEntry(i);
        std::cout << "Start Loop on event " << event << std::endl;
        int nhits = rawSvtHits->size();

        std::cout << "Number of hits: " << nhits << std::endl;
        for(int hit=0; hit < nhits; hit++){ 
            //Get hit attributes
            RawSvtHit* rawSvtHit = rawSvtHits->at(hit);
            int module = (rawSvtHit->getModule());
            int layer = (rawSvtHit->getLayer());
            std::string hwTag= mmapper_->getHwFromSw("ly"+std::to_string(layer)+"_m"+std::to_string(module));
            int channel = (int)rawSvtHit->getStrip();
            int cg = (int)channel%8;
            if(cg != select_calgroup_)
                continue;
            int svtid = svtid_map_[hwTag].at(channel); 
            std::cout << "svt id " << svtid << " channel: " << channel << std::endl;
            std::string name = hwTag+"_ch_"+std::to_string((int)channel)+"_svtid_"+std::to_string((int)svtid);

            if (pulsehistos2d_.find(name) != pulsehistos2d_.end()){
                bool exists = true;
            }
            else
                definePulseHistos(name);

            std::vector<double> adcs;
            for(int ii=0; ii < 6; ii++){
                adcs.push_back(rawSvtHit->getADCs()[ii]);
                histos2d_["adcs"]->Fill(ii,rawSvtHit->getADCs()[ii]);

                double time = 3.125*(1) + 24.0*ii;
                pulsehistos2d_[name]->Fill(time,rawSvtHit->getADCs()[ii]);
            }

            baselines_[name]->Fill(rawSvtHit->getADCs()[0]);
        }
        std::cout << "Done with event" << event << std::endl;
    }
    std::cout << "FINISHED BUILDING PULSES" << std::endl;
}

void SvtPulseFitHistos::jlab2019CalPulseScan(TTree* rawhitsTree) {

    //Establish calibration group and csel setting map based on event number
    //Calibration group that receives a pulse changes every 20,000 events, starting at calgroup = 0
    //csel setting loops from 0-8, and a "no pulse" value, with 2,000 events per setting
    std::map<std::pair<int,int>,int> cgMap = { {std::make_pair(0,20000),0},{std::make_pair(20000,40000),1},{std::make_pair(40000,60000),2},{std::make_pair(60000,80000),3},{std::make_pair(80000,100000),4},{std::make_pair(100000,120000),5},{std::make_pair(120000,140000),6},{std::make_pair(140000,160000),7}};

    std::map<std::pair<int,int>,int> cselmap = {{std::make_pair(0,2000),0}, {std::make_pair(2000,4000),1},{std::make_pair(4000,6000),2},{std::make_pair(6000,8000),3},{std::make_pair(8000,10000),4},{std::make_pair(10000,12000),5},{std::make_pair(12000,14000),6},{std::make_pair(14000,16000),7}, {std::make_pair(16000,18000),8}, {std::make_pair(18000,20000),-1} };
    double module, layer, channel, svtid, cdel, calgroup;
    double adc0, adc1, adc2, adc3, adc4, adc5;
    
    //Read rawsvthit TTree
    std::vector<RawSvtHit*>* rawSvtHits{};
    int event;
    rawhitsTree->SetBranchAddress("event", &event);
    rawhitsTree->SetBranchAddress("rawsvthits", &rawSvtHits);

    //Fit over specified calgroup only
    //by selecting corresponding events using event map
    int minevent;
    int maxevent;
    if(select_calgroup_ != -1){
        //Event calgroup and csel mapping
        std::map<std::pair<int,int>,int>::iterator it;
        int calgroup; 
        for(it=cgMap.begin(); it!=cgMap.end(); it++){
            int calgroup = it->second;
            if(calgroup == select_calgroup_){
                minevent = it->first.first;
                maxevent = (it->first.second);
                break;
            }
        }
    }

    //Loop over tree
    long nentries = rawhitsTree->GetEntries();

    //Find start and end event that matches the selected calibration group
    int start = 0;
    int end = nentries;
    if(select_calgroup_ != -1){
        start = (minevent/2)-1;
        end = maxevent/2;
    }

    int cselcount = -1;
    for(int i=start; i < end; i++){
        rawhitsTree->GetEntry(i);
        //std::cout << "Event: " << event << std::endl;

        cselcount = cselcount + 2;
        if(cselcount > 19999)
            cselcount = -1;
        if((cselcount-1)%1000==0)
            std::cout << "event " << event << std::endl;

        //Event calgroup and csel mapping
        std::map<std::pair<int,int>,int>::iterator it;
        int calgroup; 
        for(it=cgMap.begin(); it!=cgMap.end(); it++){
            int min = it->first.first;
            int max = it->first.second;
            if(event >= min && event < max){
                calgroup = it->second;
                break;
            }
        }
        int csel; 
        for(it=cselmap.begin(); it!=cselmap.end(); it++){
            int min = it->first.first;
            int max = it->first.second;
            if(cselcount >= min && cselcount < max){
                csel = it->second;
                break;
            }
        }

        //skip this group of events
        if(csel == -1)
            continue;

        int nhits = rawSvtHits->size();
        histos1d_["nhits_per_event"]->Fill(nhits);

        //Loop over all rawhits
        for(int hit=0; hit < nhits; hit++){ 
            //Get hit attributes
            RawSvtHit* rawSvtHit = rawSvtHits->at(hit);
            int module = (rawSvtHit->getModule());
            int layer = (rawSvtHit->getLayer());
            std::string hwTag= mmapper_->getHwFromSw("ly"+std::to_string(layer)+"_m"+std::to_string(module));
            int channel = (int)rawSvtHit->getStrip();
            int cg = (int)channel%8;

            //if hit channel does not belong to calgroup in this event block, skip it.
            if(cg != calgroup)
                continue;
            int svtid = svtid_map_[hwTag].at(channel); 
            std::string name = hwTag+"_ch_"+std::to_string((int)channel)+"_svtid_"+std::to_string((int)svtid);

            //This will be adjusted as histos are transformed to TGraphErrors for fitting
            if (pulsehistos2d_.find(name) != pulsehistos2d_.end()){
                bool exists = true;
            }
            else
                definePulseHistos(name);

            //Grab baseline for channel using csel=0 events
            if(csel == 0){
                baselines_[name]->Fill(rawSvtHit->getADCs()[0]);
            }
            else{
                std::vector<double> adcs;
                for(int ii=0; ii < 6; ii++){
                    adcs.push_back(rawSvtHit->getADCs()[ii]);
                    histos2d_["adcs"]->Fill(ii,rawSvtHit->getADCs()[ii]);

                    double time = getHitTime(24.0,csel,ii);
                    pulsehistos2d_[name]->Fill(time,rawSvtHit->getADCs()[ii]);
                }
            }
        }
    }
}

double SvtPulseFitHistos::getHitTime(double clockCycle, double csel, int sampleN){
    //clockCycle is in ns. 24ns for JLab
    double delay_step = clockCycle/8.0;
    double time = delay_step*(8-csel) + clockCycle*sampleN;
    return time;
}

void SvtPulseFitHistos::readPulseHistosFromFile(TFile* infile){
    std::cout << "Reading histos from file" << std::endl;
    TIter keyList(infile->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)keyList())) {
        std::string dirname;
        if((std::string)key->GetClassName() != "TDirectoryFile")
            continue;
        else {
            dirname = key->GetName();
            if(dirname.find("baselines") != std::string::npos){
                std::cout << "Found baselines dir" << std::endl;
                TDirectory* dir = (TDirectory*)infile->Get(dirname.c_str());
                TIter dirList(dir->GetListOfKeys());
                TKey *dirkey;
                while ((dirkey = (TKey*)dirList())) {
                    if((std::string)dirkey->GetClassName() != "TH1F")
                        continue;
                    std::string hname = (std::string)dirkey->GetName();
                    hname = hname.substr(0,hname.find("_baseline"));
                    std::cout << "baseline hname " << hname << std::endl;
                    TH1F* h = (TH1F*)dir->Get(dirkey->GetName());
                    //baselines_[hname] = h->GetMean();
                    //baselines_[hname] = h;
                    readbaselines_[hname] = h->GetMean();
                    delete h;
                }
            }
            else{
                TDirectory* dir = (TDirectory*)infile->Get(dirname.c_str());
                TIter dirList(dir->GetListOfKeys());
                TKey *dirkey;
                while ((dirkey = (TKey*)dirList())) {
                    if((std::string)dirkey->GetClassName() != "TH2F")
                        continue;
                    std::string hhname = (std::string)dirkey->GetName();
                    TH2F* hh = (TH2F*)dir->Get(hhname.c_str());
                    pulsehistos2d_[hhname] = hh;
                }
            }
        }
    }
    std::cout << "Finished reading " << pulsehistos2d_.size() <<" histos from file" << std::endl;
}

void SvtPulseFitHistos::buildTGraphsFromHistos(){
    //cnv25nsHistoTo24nsTGraph();
    convertHistToTGraph();
}

void SvtPulseFitHistos::fitPulses(){
    fitTGraphPulses(tGraphErrors_);
    //fit2DHistoPulses();
}

void SvtPulseFitHistos::convertHistToTGraph(){
    //Transform 2d histos to TGraph errors 
    int iter = 0;
    for (it2d it = pulsehistos2d_.begin(); it!=pulsehistos2d_.end(); ++it) {
        std::string name = it->first;

        TH2F* histo = it->second;
        int nbins = histo->GetNbinsX();
        std::vector<Double_t> means, errors, times;
        for(int i=0; i < 48; i++){
            std::string projname = name + "_bin_"+std::to_string(i+1);
            TH1F* projy = (TH1F*)histo->ProjectionY(projname.c_str(), i+1,i+1);
            double timeshift = 0;
            int entries = projy->GetEntries();
            if(entries < 1)
                continue;
            double time = (histo->GetXaxis()->GetBinCenter(i+1)) - timeshift;
            double mean = projy->GetMean();
            double err = projy->GetStdDev();
            means.push_back(mean);
            errors.push_back(err);
            times.push_back(time);

            delete projy;
        }

        TGraphErrors* tgr = new TGraphErrors((int)means.size(),&(times[0]),&(means[0]), 0, &(errors[0]) ); 
        tgr->SetName((name+"_tgr").c_str());
        tgr->SetTitle((name+";time (ns);ADC").c_str());
        tGraphErrors_[name] = tgr;
        iter++;
    }
}

void SvtPulseFitHistos::checkPulseQuality(TGraphErrors* tgraph, bool &badPulse){

    if(year_ == 2019){
        int npoints = tgraph->GetN();
        int badErr = 0;
        for(int i=0; i < npoints; i++){
            double error = tgraph->GetErrorY(i);
            if(error >=350){
                badErr++;
            }
        }

        if(badErr > 5)
            badPulse = true;
        
        if(tgraph->GetPointY(11) < 500 + tgraph->GetPointY(0))
            badPulse = true;
    }
    else{
        int npoints = tgraph->GetN();
        for(int i=1; i < npoints-1; i++){
            if(tgraph->GetPointY(i) < tgraph->GetPointY(i-1) && tgraph->GetPointY(i) < tgraph->GetPointY(i+1)){
               badPulse = true; 
            }
        }
        if(tgraph->GetPointY(3) < 500.0 + tgraph->GetPointY(0))
            badPulse = true;
        return;
    }
}

void SvtPulseFitHistos::fit2DHistoPulses(){
    int iiter = 0;
    for (it2d it = pulsehistos2d_.begin(); it!=pulsehistos2d_.end(); ++it) {
        if(iiter>10)
            break;

        TH2F* pulsehisto2d = it->second;

        //Get channel names and info
        //std::string s = pulsehisto2d->GetName();
        std::string s = it->first;
        std::cout << "Fitting 2dhisto : " << s << std::endl;
        int svtid = std::stoi(s.substr(s.find("svtid_")+6,s.size()-1));
        std::string hwTag = s.substr(0,s.find("_"));
        int channel = std::stoi(s.substr(s.find("ch_")+3,s.find("_s")-5));
        std::string sw = mmapper_->getSwFromHw(hwTag);
        int layer = std::stoi(sw.substr(sw.find("ly")+2,sw.find("_")-1));
        int module = std::stoi(sw.substr(sw.find("m")+1,sw.size()-1));

        //amplitude seed
        double maxamp=0;
        for (int i = 0; i < pulsehisto2d->GetNbinsX(); i++) {
            TH1F* proj = (TH1F*)pulsehisto2d->ProjectionY(((std::string)pulsehisto2d->GetName()+"_proj").c_str(), i+1, i+1);
            double a = proj->GetMean();
            if(a > maxamp){
                maxamp = a;
            }
            delete proj;
        }

        //fix baseline using csel = 0 events
        double baseline = baselines_[s]->GetMean();

        //APV25 Channel Response Fit Function
        TF1* fitfunc = fourPoleFitFunction();
        //Random seed generator
        TRandom *r1 = new TRandom();

        //Seed fit parameters
        double bestchi2ndf = 99999999.0;
        double t0;
        double tau1;
        double tau2;
        double amp;
        bool nanfit = false;
        int maxiter = 40;
        int iter = 0;
        while(iter < maxiter){
            double rtau1;
            double rtau2;
            double rt0;
            double ramp; 

            //Slim sensor random fit param seeds
            if (svtid < 4096){
                rtau1 = r1->Gaus(53.,2.);
                rtau2 = r1->Gaus(9.,2.);
                rt0 = r1->Uniform(-10.,10.);
                ramp = r1->Gaus(maxamp,2.);
            }
            //Non-slim sensor random fit param seeds
            else {
                rtau1 = r1->Gaus(55.,2.);
                rtau2 = r1->Gaus(12.,2.);
                rt0 = r1->Uniform(-10.,10.);
                ramp = r1->Gaus(maxamp,2);
            }

            //set function seeds
            fitfunc->SetParameter(1,rtau1);
            fitfunc->SetParameter(2,rtau2);
            fitfunc->SetParameter(0,rt0);
            fitfunc->SetParameter(3,ramp);
            fitfunc->FixParameter(4,baseline);

            //Fit using random seeds
            pulsehisto2d->Fit(fitfunc, "q");

            //Get fit quality via chi2
            double chi2 = fitfunc->GetChisquare();
            double ndf = (double) fitfunc->GetNDF();
            if(ndf < 1){
                iter++;
                continue;
            }
            double chi2ndf = chi2/ndf;

            //Keep fit with best chi2
            if(chi2ndf < bestchi2ndf){
                bestchi2ndf = chi2ndf;
                tau1 = rtau1;
                tau2 = rtau2;
                t0 = rt0;
                amp = ramp;
                ndf = ndf;
            }

            iter++;

            //If bestchi2ndf < 1, fit should be good enough to stop iteration
            if(bestchi2ndf < 1)
                break;
        }

        //Refit with parameters set to values derived from fit with best chi2
        fitfunc->SetParameter(1,tau1);
        fitfunc->SetParameter(2,tau2);
        fitfunc->SetParameter(0,t0);
        fitfunc->SetParameter(3,amp);
        fitfunc->FixParameter(4,baseline);
        pulsehisto2d->Fit(fitfunc, "q");

        //Get final fit parameters
        t0 = fitfunc->GetParameter(0);
        tau1 = fitfunc->GetParameter(1);
        tau2 = fitfunc->GetParameter(2);
        amp = fitfunc->GetParameter(3);
        baseline = fitfunc->GetParameter(4);
        double chi2 = fitfunc->GetChisquare();
        double ndf = fitfunc->GetNDF();
        if(ndf < 1){
            nanfit = true;
        }

        double t0err = fitfunc->GetParError(0);
        double tau1err = fitfunc->GetParError(1);
        double tau2err = fitfunc->GetParError(2);
        double amperr = fitfunc->GetParError(3);

        pulsehisto2d->Draw();
        iiter++;
    }
}

void SvtPulseFitHistos::fitTGraphPulses(std::map<std::string,TGraphErrors*> tgraphs){
    std::cout << "FITTING TGRAPHS" << std::endl;

    for (it_tgr it = tgraphs.begin(); it!=tgraphs.end(); ++it) {

        TGraphErrors* tgraph = it->second;
        //Get channel names and info
        //std::string name = tgraph->GetName();
        std::string name = it->first;
        std::cout << "Fitting tgraph : " << name << std::endl;
        int svtid = std::stoi(name.substr(name.find("svtid_")+6,name.size()-1));
        std::string hwTag = name.substr(0,name.find("_"));
        int channel = std::stoi(name.substr(name.find("ch_")+3,name.find("_s")-5));
        std::string sw = mmapper_->getSwFromHw(hwTag);
        int layer = std::stoi(sw.substr(sw.find("ly")+2,sw.find("_")-1));
        int module = std::stoi(sw.substr(sw.find("m")+1,sw.size()-1));

        //amplitude seed
        double maxamp=0;
        std::cout << pulsehistos2d_[name]->GetName() << std::endl;
        for (int i = 0; i < pulsehistos2d_[name]->GetNbinsX(); i++) {
            TH1F* proj = (TH1F*)pulsehistos2d_[name]->ProjectionY(((std::string)pulsehistos2d_[name]->GetName()+"_proj").c_str(), i+1, i+1);
            double a = proj->GetMean();
            if(a > maxamp){
                maxamp = a;
            }
            delete proj;
        }

        //double baseline = baselines_[name]->GetMean();
        double baseline = readbaselines_[name];

        //APV25 Channel Response Fit Function
        TF1* fitfunc = fourPoleFitFunction();
        //Random seed generator
        TRandom *r1 = new TRandom();

        //fit quality
        bool nanfit = false;
        bool badPulse = false;

        //refit nanfits
        int notnancount = 0;
        int maxnaniter = 5000;

        //Best fit parameters
        double bestchi2ndf = 99999999.0;
        double bestchi2;
        double bestndf;
        double bestt0;
        double besttau1;
        double besttau2;
        double bestamp;
        double bestt0err;
        double besttau1err;
        double besttau2err;
        double bestamperr;

        //iter
        int maxiter = 200;
        int iter = 0;
        while( (iter < maxiter) || (notnancount < 5 && iter < maxnaniter) ){
            //random generated fit seeds
            double rtau1;
            double rtau2;
            double rt0;
            double ramp; 

            //Fit result params
            double chi2ndf;
            double chi2;
            double ndf;
            double t0;
            double tau1;
            double tau2;
            double amp;
            double t0err;
            double tau1err;
            double tau2err;
            double amperr;

            std::cout << "Fit iteration " << iter << std::endl;

            //Slim sensor random fit param seeds
            if (svtid < 4096){
                rtau1 = r1->Gaus(53.,1.);
                rtau2 = r1->Gaus(8.,1.);
                //rt0 = r1->Uniform(-10.,-1.);
                rt0 = r1->Gaus(20.,1);
                ramp = r1->Gaus(maxamp,1);
            }
            //Non-slim sensor random fit param seeds
            else {
                rtau1 = r1->Gaus(53.,1);
                rtau2 = r1->Gaus(11.,1);
                //rt0 = r1->Uniform(-10.,-1.);
                rt0 = r1->Gaus(20.,1);
                ramp = r1->Gaus(maxamp,1);
            }

            //set function seeds
            fitfunc->SetParameter(1,rtau1);
            fitfunc->SetParameter(2,rtau2);
            fitfunc->SetParameter(0,rt0);
            fitfunc->SetParameter(3,ramp);
            fitfunc->FixParameter(4,baseline);

            //Fit using random seeds
            tgraph->Fit(fitfunc, "q");

            //Get fit result parameters
            t0 = fitfunc->GetParameter(0);
            tau1 = fitfunc->GetParameter(1);
            tau2 = fitfunc->GetParameter(2);
            amp = fitfunc->GetParameter(3);
            baseline = fitfunc->GetParameter(4);
            chi2 = fitfunc->GetChisquare();
            ndf = (double)fitfunc->GetNDF();
            t0err = fitfunc->GetParError(0);
            tau1err = fitfunc->GetParError(1);
            tau2err = fitfunc->GetParError(2);
            amperr = fitfunc->GetParError(3);

            std::cout << "Chi2 of iteration " << iter << ": " << chi2 << std::endl;
            std::cout << "NDF of iteration "<< iter << ": " << ndf << std::endl;

            //Check if fit result is NAN
            if (TMath::IsNaN(t0err) || TMath::IsNaN(tau1err) || TMath::IsNaN(tau2err) || TMath::IsNaN(amperr)){
                nanfit = true;
                std::cout << "Found NAN fit in iter" << iter << std::endl;
                iter++;
                continue;
            }
            else{
               std::cout << "No NAN in this fit iteration" << iter << std::endl;
               nanfit = false;
               notnancount++; 
               std::cout << "Not NAN count = " << notnancount << std::endl;
            }
            if(ndf < 1){
                std::cout << "NDF < 1! BAD FIT!" << std::endl;
                iter++;
                continue;
            }
            //Get chi2 of this fit iteration
            chi2ndf = chi2/ndf;

            std::cout << "Best chi2 so far is " << bestchi2ndf << std::endl;
            //Keep fit with best chi2
            if(chi2ndf < bestchi2ndf){
                bestchi2ndf = chi2ndf;
                std::cout << "Best chi2 updated to " << bestchi2ndf << " on iteration " << iter << std::endl;
                besttau1 = tau1;
                besttau2 = tau2;
                bestt0 = t0;
                bestamp = amp;
                besttau1err = tau1err;
                besttau2err = tau2err;
                bestt0err = t0err;
                bestamperr = amperr;
                bestndf = ndf;
                bestchi2 = chi2;
            }
            iter++;
        }

        std::cout << "svtid " << svtid << "after iteration" <<  std::endl;
        std::cout << "t0err: " << bestt0err << " | tau1err: " << besttau1err << " | tau2err: " << besttau2err << " | amperr: " << bestamperr << std::endl;

        //check final fit for NAN results
        if (TMath::IsNaN(bestt0err) || TMath::IsNaN(besttau1err) || TMath::IsNaN(besttau2err) || TMath::IsNaN(bestamperr)){
            nanfit = true;
            std::cout << "Final fit has NAN errors" << std::endl;
            bestt0err = -99999;
            besttau1err = -99999;
            besttau2err = -99999;
            bestamperr = -99999;
        }
        else{
            nanfit = false;
            std::cout << "Final Fit has NO NAN ERRORS" << std::endl;
        }

        //Check pulse to make sure its not a "bad pulse"
        checkPulseQuality(tgraph, badPulse);

        std::cout << "svtid " << svtid << " setting tuple" <<  std::endl;
        std::cout << "t0err: " << bestt0err << " | tau1err: " << besttau1err << " | tau2err: " << besttau2err << " | amperr: " << bestamperr << std::endl;
        //add fit parameter results to tuple for export
        rawhitfits_tup_->setVariableValue("hwTag", hwTag);
        rawhitfits_tup_->setVariableValue("channel", channel);
        rawhitfits_tup_->setVariableValue("svtid", svtid);
        rawhitfits_tup_->setVariableValue("module", module);
        rawhitfits_tup_->setVariableValue("layer", layer);
        rawhitfits_tup_->setVariableValue("nanfit", nanfit);
        rawhitfits_tup_->setVariableValue("badPulse", badPulse);
        rawhitfits_tup_->setVariableValue("t0", bestt0);
        rawhitfits_tup_->setVariableValue("tau1", besttau1);
        rawhitfits_tup_->setVariableValue("tau2", besttau2);
        rawhitfits_tup_->setVariableValue("amp", bestamp);
        rawhitfits_tup_->setVariableValue("baseline", baseline);
        rawhitfits_tup_->setVariableValue("chi2", bestchi2);
        rawhitfits_tup_->setVariableValue("ndf", bestndf);
        rawhitfits_tup_->setVariableValue("t0err", bestt0err);
        rawhitfits_tup_->setVariableValue("tau1err", besttau1err);
        rawhitfits_tup_->setVariableValue("tau2err", besttau2err);
        rawhitfits_tup_->setVariableValue("amperr", bestamperr);
        rawhitfits_tup_->setVariableValue("integralNorm", SvtPulseFitHistos::getAmplitudeIntegralNorm(besttau1, besttau2));

        rawhitfits_tup_->fill();

        //Refit with the best parameters to save TGraph Fit
        fitfunc->FixParameter(1,besttau1);
        fitfunc->FixParameter(2,besttau2);
        fitfunc->FixParameter(0,bestt0);
        fitfunc->FixParameter(3,bestamp);
        fitfunc->FixParameter(4,baseline);
        tgraph->Fit(fitfunc, "q");

        double t0err = fitfunc->GetParError(0);
        double tau1err = fitfunc->GetParError(1);
        double tau2err = fitfunc->GetParError(2);
        double amperr = fitfunc->GetParError(3);

        std::cout << "svtid " << svtid << "post tuple fit" <<  std::endl;
        std::cout << "t0err: " << t0err << " | tau1err: " << tau1err << " | tau2err: " << tau2err << " | amperr: " << amperr << std::endl;

        //Fill Histograms
        if(bestndf > 0){
            histos1d_["chi2"]->Fill(bestchi2ndf);
        }
        histos1d_["tau1_v_id"]->SetBinContent(svtid+1,besttau1);
        histos1d_["tau2_v_id"]->SetBinContent(svtid+1,besttau2);
        histos2d_["t0_v_amp"]->Fill(fitfunc->GetParameter(0),bestamp,1.0);
        histos2d_["tau1_v_tau2"]->Fill(besttau1,besttau2,1.0);

        tgraph->Draw();
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

void SvtPulseFitHistos::saveHistos(TFile* outFile){
    outFile->cd();

    //Save 2d pulse histos
    std::cout << "SAVING PULSE 2D HISTOS: " << pulsehistos2d_.size() << std::endl;
    for (it2d it = pulsehistos2d_.begin(); it!=pulsehistos2d_.end(); ++it) {

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
        if(!exists){
            gDirectory->mkdir(dirname.c_str());
        }

        outFile->cd(dirname.c_str());
        it->second->Write();
    }

    //Save TGraphErrors
    std::cout << "SAVING TGRAPH ERRORS: " << tGraphErrors_.size() << std::endl;
    for (it_tgr it = tGraphErrors_.begin(); it!=tGraphErrors_.end(); ++it) {

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
        if(!exists){
            gDirectory->mkdir(dirname.c_str());
        }

        outFile->cd(dirname.c_str());
        it->second->Write();
    }

    //save baseline histos
    outFile->cd();
    gDirectory->mkdir("baselines");
    outFile->cd("baselines");
    for(it1d it = baselines_.begin(); it!=baselines_.end(); it++){
        it->second->Write();
    }

    //save 1d histos
    std::cout << "SAVING 1D HISTOS" << histos1d_.size() << std::endl;
    outFile->cd();
    for (it1d it = histos1d_.begin(); it!=histos1d_.end(); it++) {
        it->second->Write();
    }
    //save 2d histos
    std::cout << "SAVING 2D HISTOS" << histos2d_.size() << std::endl;
    outFile->cd();
    typedef std::map<std::string, TH2F*>::iterator it2d;
    for (it2d it = histos2d_.begin(); it!=histos2d_.end(); it++) {
        it->second->Write();
    }
}
