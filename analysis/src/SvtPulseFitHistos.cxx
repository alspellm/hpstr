#include "SvtPulseFitHistos.h"
#include <math.h>
#include "TCanvas.h"
#include "TRandom.h"

SvtPulseFitHistos::SvtPulseFitHistos(const std::string& inputName, ModuleMapper* mmapper) {
    m_name = inputName;
    mmapper_ = mmapper;
    svtid_map_ = mmapper_->buildChannelSvtIDMap();
}

SvtPulseFitHistos::~SvtPulseFitHistos() {
}

void SvtPulseFitHistos::initHistos(){
    //1d histos
    TH1F* hpe_h = new TH1F("nhits_per_event","nhits_per_event",25000,0,25000);
    TH1F* chi2_h = new TH1F("bestchi2","best_chi2;chi2;entries",1000,0,10);
    TH1F* tau1_v_id = new TH1F("tau1_v_id","tau1_v_id",25000,-0.5,24999.5);
    TH1F* tau2_v_id = new TH1F("tau2_v_id","tau2_v_id",25000,-0.5,24999.5);

    histos1d_["nhits_per_event"] = hpe_h;
    histos1d_["chi2"] = chi2_h;
    histos1d_["tau1_v_id"] = tau1_v_id;
    histos1d_["tau2_v_id"] = tau2_v_id;

    //2d histos
    TH2F* adcs_hh = new TH2F("adcs","adcs;sample N;ADC",6,0,6,5000,0,10000);
    TH2F* t0_v_amp_hh = new TH2F("t0_v_amp","t0_v_amp;t0 (ns);Amp (adc)",100,0,100,4000,0,400000);
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
    //Pulse Histos are binned under the assumption of a 25ns clock (true for UCSC testboard data)
    //If using JLab DAQ pulse data, the clock is 24ns, and the time bins need to be adjusted...this is done
    //by the method SvtPulseFitHistos::adjustClock25nsTo24ns
    TH2F* prof = new TH2F(name.c_str(),name.c_str(), 49,-1.5625,151.5625,10000,0,10000);
    pulsehistos2d_[name] = prof;

    //store channel baselines in histogram
    TH1F* baseline = new TH1F((name+"_baseline").c_str(),name.c_str(),10000,0,10000);
    baselines_[name] = baseline;
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
    rawhitsTree->Print();
    std::vector<RawSvtHit*>* rawSvtHits{};
    int event;
    rawhitsTree->SetBranchAddress("event", &event);
    rawhitsTree->SetBranchAddress("rawsvthits", &rawSvtHits);

    //Loop over tree
    long nentries = rawhitsTree->GetEntries();
    std::cout << "TTree N Entries: " << nentries << std::endl;
    int cselcount = -1;
    for(long i=0; i < nentries; i++){
        rawhitsTree->GetEntry(i);
        if((cselcount-1)%1000==0)
            std::cout << "event " << event << std::endl;
        cselcount = cselcount + 2;
        if(cselcount > 19999)
            cselcount = -1;
        /*
        //Cut on event
        if(event > 20000)
            break;
            */
        
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
        if(csel == -1)
            continue;

        //std::cout << "calgroup = " << calgroup << std::endl;
        //std::cout << "csel = " << csel << std::endl;
        //
        int nhits = rawSvtHits->size();
        //std::cout << "nhits = " << nhits << std::endl;
        histos1d_["nhits_per_event"]->Fill(nhits);

        for(int hit=0; hit < nhits; hit++){ 
            //Get hit attributes
            RawSvtHit* rawSvtHit = rawSvtHits->at(hit);
            int module = (rawSvtHit->getModule());
            int layer = (rawSvtHit->getLayer());
            std::string hwTag= mmapper_->getHwFromSw("ly"+std::to_string(layer)+"_m"+std::to_string(module));
            int channel = (int)rawSvtHit->getStrip();
            int cg = (int)channel%8;
            if(cg != calgroup)
                continue;
            int svtid = svtid_map_[hwTag].at(channel); 
            std::string name = hwTag+"_ch_"+std::to_string((int)channel)+"_svtid_"+std::to_string((int)svtid);

            /*
            //dev cut on channels
            if(svtid > 10)
                continue;
                */

            //Channel histograms defined below use incorrect clock binning (25.ns, should be 24.ns)
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

                    double time = 3.125*(8-csel) + 25.0*ii;
                    pulsehistos2d_[name]->Fill(time,rawSvtHit->getADCs()[ii]);
                }
            }
        }
    }

    //Take UCSC Testbaord 25ns 2d histos and translate them into JLab DAQ 24ns TGraphErrors 
    adjustClock25nsTo24ns();
}

void SvtPulseFitHistos::adjustClock25nsTo24ns(){
    //Transform 2d histos to TGraph errors with clock adjusted from 25ns to 24 ns
    for (it2d it = pulsehistos2d_.begin(); it!=pulsehistos2d_.end(); ++it) {
        std::string name = it->first;
        TH2F* histo = it->second;
        int nbins = histo->GetNbinsX();
        //std::cout << "nbins = " << nbins << std::endl;
        std::vector<Double_t> means, errors, times;
        //std::vector<double> errors;
        //std::vector<double> hhtimes;
        for(int i=0; i < 48; i++){
            std::string projname = name + "_bin_"+std::to_string(i+1);
            TH1F* projy = (TH1F*)histo->ProjectionY(projname.c_str(), i+1,i+1);
            double timeshift = floor((i+1)/9);
            int entries = projy->GetEntries();
            if(entries < 1)
                continue;
            //std::cout << "time before shift = " << (histo->GetXaxis()->GetBinCenter(i+1)) << std::endl;
            //std::cout << "timeshift = " << timeshift << std::endl;
            double time = (histo->GetXaxis()->GetBinCenter(i+1)) - timeshift;
            //std::cout << "time after shift = " << time << std::endl;
            double mean = projy->GetMean();
            double err = projy->GetStdDev();
            means.push_back(mean);
            errors.push_back(err);
            times.push_back(time);

            delete projy;
        }
        /*
        double adjustedTimes[(int)hhtimes.size()];

        for(int i=0; i<8; i++){
            int iter = 0;
            for(int ii=i; ii<48; ii+8){
                double t_adjusted = hhtimes.at(ii) - iter; 
                adjustedTimes[ii] = t_adjusted;
                iter = iter + 1;
            }
        }
        */
        TGraphErrors* tgr = new TGraphErrors((int)means.size(),&(times[0]),&(means[0]), 0, &(errors[0]) ); 
        tgr->SetName((name+"_tgr").c_str());
        tgr->SetTitle((name+";time (ns);ADC").c_str());
        tgrapherrs_[name] = tgr;
    }

    //Fit TGraphErrors
    fitTGraphPulses();
    //Fit histos
    fit2DHistoPulses();
}

void SvtPulseFitHistos::fit2DHistoPulses(){
    for (it2d it = pulsehistos2d_.begin(); it!=pulsehistos2d_.end(); ++it) {

        TH2F* pulsehisto2d = it->second;

        //Get channel names and info
        //std::string s = pulsehisto2d->GetName();
        std::string s = it->first;
        int svtid = std::stoi(s.substr(s.find("svtid_")+6,s.size()-1));
        std::string hwTag = s.substr(0,s.find("_"));
        int channel = std::stoi(s.substr(s.find("ch_")+3,s.find("_s")-5));
        std::string sw = mmapper_->getSwFromHw(hwTag);
        int layer = std::stoi(sw.substr(sw.find("ly")+2,sw.find("_")-1));
        int module = std::stoi(sw.substr(sw.find("m")+1,sw.size()-1));

        //amplitude seed
        //TH1F* projy = (TH1F*)pulsehisto2d->ProjectionY((s+"tmp_projy").c_str(), 11,11);
        double maxamp;
        for (int i = 0; i < pulsehisto2d->GetNbinsX(); i++) {
            double a = pulsehisto2d->GetBinContent(i+1);
            if(a > maxamp){
                maxamp = a;
            }
        }

        //delete projy;
        //fix baseline using csel = 0 events
        double baseline = baselines_[s]->GetMean();

        //APV25 Channel Response Fit Function
        TF1* fitfunc = fourPoleFitFunction();
        //Random seed generator
        TRandom *r1 = new TRandom();

        //Seed fit parameters
        double bestchi2 = 99999999.0;
        double t0;
        double tau1;
        double tau2;
        double amp;
        bool nanfit = false;
        int maxiter = 30;
        int iter = 0;
        while(iter < maxiter){
            double rtau1;
            double rtau2;
            double rt0;
            double ramp; 

            //Slim sensor random fit param seeds
            if (svtid < 4096){
                rtau1 = r1->Gaus(53.,3.);
                rtau2 = r1->Gaus(9.,3.);
                rt0 = r1->Uniform(0.,30.);
                ramp = r1->Gaus(maxamp,3.);
            }
            //Non-slim sensor random fit param seeds
            else {
                rtau1 = r1->Gaus(55.,3.);
                rtau2 = r1->Gaus(12.,3.);
                rt0 = r1->Uniform(0.,30.);
                ramp = r1->Gaus(maxamp,3);
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
            if(chi2ndf < bestchi2){
                bestchi2 = chi2ndf;
                tau1 = rtau1;
                tau2 = rtau2;
                t0 = rt0;
                amp = ramp;
                ndf = ndf;
            }

            iter++;

            //If bestchi2 < 1, fit should be good enough to stop iteration
            if(bestchi2 < 1)
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
    }

}

void SvtPulseFitHistos::fitTGraphPulses(){

    for (it_tgr it = tgrapherrs_.begin(); it!=tgrapherrs_.end(); ++it) {

        TGraphErrors* tgraph = it->second;

        //Get channel names and info
        //std::string s = tgraph->GetName();
        std::string s = it->first;
        int svtid = std::stoi(s.substr(s.find("svtid_")+6,s.size()-1));
        std::string hwTag = s.substr(0,s.find("_"));
        int channel = std::stoi(s.substr(s.find("ch_")+3,s.find("_s")-5));
        std::string sw = mmapper_->getSwFromHw(hwTag);
        int layer = std::stoi(sw.substr(sw.find("ly")+2,sw.find("_")-1));
        int module = std::stoi(sw.substr(sw.find("m")+1,sw.size()-1));

        //amplitude seed
        //TH1F* projy = (TH1F*)pulsehistos2d_[s]->ProjectionY((s+"tmp_projy").c_str(), 11,11);
        //double maxamp = projy->GetMean();
        double maxamp;
        for (int i = 0; i < pulsehistos2d_[s]->GetNbinsX(); i++) {
            double a = pulsehistos2d_[s]->GetBinContent(i+1);
            if(a > maxamp){
                maxamp = a;
            }
        }

        //double maxamp = tgraph->GetMaximum();
        std::cout << "TGRAPH GET MAXIMUM: " << maxamp << std::endl;
        //fix baseline using csel = 0 events
        double baseline = baselines_[s]->GetMean();

        //APV25 Channel Response Fit Function
        TF1* fitfunc = fourPoleFitFunction();
        //Random seed generator
        TRandom *r1 = new TRandom();

        //Seed fit parameters
        double bestchi2 = 99999999.0;
        double t0;
        double tau1;
        double tau2;
        double amp;
        bool nanfit = false;
        int maxiter = 30;
        int iter = 0;
        while(iter < maxiter){
            double rtau1;
            double rtau2;
            double rt0;
            double ramp; 

            //Slim sensor random fit param seeds
            if (svtid < 4096){
                rtau1 = r1->Gaus(53.,3.);
                rtau2 = r1->Gaus(9.,3.);
                rt0 = r1->Uniform(0.,30.);
                ramp = r1->Gaus(maxamp,3.);
            }
            //Non-slim sensor random fit param seeds
            else {
                rtau1 = r1->Gaus(55.,3.);
                rtau2 = r1->Gaus(12.,3.);
                rt0 = r1->Uniform(0.,30.);
                ramp = r1->Gaus(maxamp,3);
            }

            //set function seeds
            fitfunc->SetParameter(1,rtau1);
            fitfunc->SetParameter(2,rtau2);
            fitfunc->SetParameter(0,rt0);
            fitfunc->SetParameter(3,ramp);
            fitfunc->FixParameter(4,baseline);

            //Fit using random seeds
            tgraph->Fit(fitfunc, "q");

            //Get fit quality via chi2
            double chi2 = fitfunc->GetChisquare();
            double ndf = (double) fitfunc->GetNDF();
            if(ndf < 1){
                iter++;
                continue;
            }
            double chi2ndf = chi2/ndf;

            //Keep fit with best chi2
            if(chi2ndf < bestchi2){
                bestchi2 = chi2ndf;
                tau1 = rtau1;
                tau2 = rtau2;
                t0 = rt0;
                amp = ramp;
                ndf = ndf;
            }

            iter++;

            //If bestchi2 < 1, fit should be good enough to stop iteration
            if(bestchi2 < 1)
                break;
        }

        //Refit with parameters set to values derived from fit with best chi2
        fitfunc->SetParameter(1,tau1);
        fitfunc->SetParameter(2,tau2);
        fitfunc->SetParameter(0,t0);
        fitfunc->SetParameter(3,amp);
        fitfunc->FixParameter(4,baseline);
        tgraph->Fit(fitfunc, "q");

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

        if (TMath::IsNaN(fitfunc->GetParError(0)) || TMath::IsNaN(fitfunc->GetParError(1)) || TMath::IsNaN(fitfunc->GetParError(2)) || TMath::IsNaN(fitfunc->GetParError(3))){
            nan_channels_++;
            nanfit = true;
        }

        rawhitfits_tup_->setVariableValue("hwTag", hwTag);
        rawhitfits_tup_->setVariableValue("channel", channel);
        rawhitfits_tup_->setVariableValue("svtid", svtid);
        rawhitfits_tup_->setVariableValue("module", module);
        rawhitfits_tup_->setVariableValue("layer", layer);
        rawhitfits_tup_->setVariableValue("nanfit", nanfit);
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

        if(nanfit == false){
            //Fill Histograms
            histos1d_["chi2"]->Fill(fitfunc->GetChisquare()/fitfunc->GetNDF());
            histos1d_["tau1_v_id"]->SetBinContent(svtid+1,fitfunc->GetParameter(1));
            histos1d_["tau2_v_id"]->SetBinContent(svtid+1,fitfunc->GetParameter(2));

            histos2d_["t0_v_amp"]->Fill(fitfunc->GetParameter(0),fitfunc->GetParameter(3),1.0);
            histos2d_["tau1_v_tau2"]->Fill(fitfunc->GetParameter(1),fitfunc->GetParameter(2),1.0);
        }

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
    std::cout << "SAVING TGRAPH ERRORS: " << tgrapherrs_.size() << std::endl;
    for (it_tgr it = tgrapherrs_.begin(); it!=tgrapherrs_.end(); ++it) {

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
    std::cout << "Number of NAN channel fits: " << nan_channels_ << std::endl;
}

