#include "SvtPulseFitHistos.h"
#include <math.h>
#include "TCanvas.h"

SvtPulseFitHistos::SvtPulseFitHistos(const std::string& inputName, ModuleMapper* mmapper) {
    m_name = inputName;
    mmapper_ = mmapper;
    svtid_map_ = mmapper_->buildChannelSvtIDMap();
}

SvtPulseFitHistos::~SvtPulseFitHistos() {
}

double SvtPulseFitHistos::GetHitTime(int sample_number, int cdel) {
    double time = 24.0*sample_number;
    return time;
}

void SvtPulseFitHistos::buildRawSvtHitsTuple(std::vector<RawSvtHit*> *rawSvtHits_, FlatTupleMaker* rawhits_tup_) {

    bool debug = false;
    int nhits = rawSvtHits_->size();
    std::cout << "NHITS = " << nhits << std::endl;
    for (int i = 0; i < nhits; i++)
    {
        RawSvtHit* rawSvtHit = rawSvtHits_->at(i);
        int module = (rawSvtHit->getModule());
        int layer = (rawSvtHit->getLayer());
        std::string hwTag= mmapper_->getHwFromSw("ly"+std::to_string(layer)+"_m"+std::to_string(module));
        float channel = rawSvtHit->getStrip();
        int svtid = mmapper_->getSvtIDFromHWChannel(channel, hwTag, svtid_map_); 
        int calgroup = (int)channel%8;
        //Hard coded for now. Need to find way to determine this from data
        int cdel = 1;

        if(debug){
            std::cout << "module: " << module << " layer: " << layer << " channel: " << channel << 
                " svtid: " << svtid << " calgroup: " << calgroup << std::endl;
        }

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
}

void SvtPulseFitHistos::defineTProfiles(int maxchannels){
    for(int i = 0; i < maxchannels; i++){
        std::string name = "svtid_"+std::to_string(i);
        TProfile* prof = new TProfile(name.c_str(),name.c_str(), 73,-1,145,2000,8000);
        tprofiles_[name] = prof;
    }
}

void SvtPulseFitHistos::defineTProfile(std::string name) {
    TProfile* prof = new TProfile(name.c_str(),name.c_str(), 73,-1,145,2000,8000);
    tprofiles_[name] = prof;
}

void SvtPulseFitHistos::fitRawHitPulses(TTree* rawhittree) {

    //define TProfiles for all channels
    //defineTProfiles(200);


    double module, layer, channel, svtid, cdel, calgroup;
    double adc0, adc1, adc2, adc3, adc4, adc5;
    std::string hwTag;

    rawhittree->SetBranchAddress("svtid", &svtid);
    rawhittree->SetBranchAddress("adc0", &adc0);
    rawhittree->SetBranchAddress("adc1", &adc1);
    rawhittree->SetBranchAddress("adc2", &adc2);
    rawhittree->SetBranchAddress("adc3", &adc3);
    rawhittree->SetBranchAddress("adc4", &adc4);
    rawhittree->SetBranchAddress("adc5", &adc5);
    rawhittree->SetBranchAddress("cdel", &cdel);

    /*
    for (ittp it = tprofiles_.begin(); it!=tprofiles_.end(); ++it) {
        std::cout << "TP NAMES: " <<  it->first << std::endl;
        //it->second->Fill(100.0,100.0,1.0);
        tprofiles_[it->first]->Fill(100.0,100.0,1.0);
    }
    */

    long nentries = rawhittree->GetEntries();
    for(long i=0; i < nentries; i++){
        rawhittree->GetEntry(i);
        //Build adcs vector
        std::vector<double> adcs;
        adcs.push_back(adc0);
        adcs.push_back(adc1);
        adcs.push_back(adc2);
        adcs.push_back(adc3);
        adcs.push_back(adc4);
        adcs.push_back(adc5);

        std::string name = "svtid_"+std::to_string((int)svtid);

        //check if tprofile exists
        if (tprofiles_.find(name) != tprofiles_.end()) {
            bool exists = true;
        }
        else{
           defineTProfile(name); 
        }

        for(int t=0; t < 6; t++){
            //std::cout << "adc" << t << ": " << adc << " at time " << time << std::endl;
            double adc = adcs.at(t);
            double time = GetHitTime(t, cdel);
            tprofiles_[name]->Fill(time,adc,1.0);
        }
    }

    for (ittp it = tprofiles_.begin(); it!=tprofiles_.end(); ++it) {
        fitPulse(it->second);
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

void SvtPulseFitHistos::fitPulse(TProfile* tprofile){

    double t0 = 0.0;
    double tau1 = 55.0;
    double tau2 = 8.0;
    double amp = 0.0;
    double baseline = 0.0;

    //Get amplitude and baseline seed
    double maxamp = 0.0;
    double count = 0.0;
    for (int i = 0; i < tprofile->GetNbinsX(); i++) {
        std::cout << tprofile->GetName() << std::endl;
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

    baseline = baseline/count;
    amp = maxamp;

    //std::cout << "baseline seed: " << baseline << std::endl;
    //std::cout << "amp seed: " << amp << std::endl;

    TF1* fitfunc = fourPoleFitFunction();
    fitfunc->SetParameter(0,t0);
    fitfunc->SetParameter(1,tau1);
    fitfunc->SetParameter(2,tau2);
    fitfunc->SetParameter(3,amp);
    fitfunc->FixParameter(4,baseline);
    
    tprofile->Fit(fitfunc, "q");
    tprofile->Draw();
    tprofile->Write();

}

void SvtPulseFitHistos::saveTProfiles(TFile* outF, std::string folder) {
    outF->cd();
    for (ittp it = tprofiles_.begin(); it!=tprofiles_.end(); ++it) {

        it->second->Write();
    }
}

void SvtPulseFitHistos::FillHistogramsByHw(std::vector<RawSvtHit*> *rawSvtHits_,float weight) {

    int nhits = rawSvtHits_->size();
    std::vector<std::string> hybridStrings={};
    std::string histokey;
    if(Event_number%10000 == 0) std::cout << "Event: " << Event_number 
        << " Number of RawSvtHits: " << nhits << std::endl;

    //Populates histograms for each hybrid channel
    for (int i = 0; i < nhits; i++)
    {
        RawSvtHit* rawSvtHit = rawSvtHits_->at(i);
        auto mod = std::to_string(rawSvtHit->getModule());
        auto lay = std::to_string(rawSvtHit->getLayer());
        std::string hwTag= mmapper_->getHwFromSw("ly"+lay+"_m"+mod);
        float channel = rawSvtHit->getStrip();
        int svtid = mmapper_->getSvtIDFromHWChannel(channel, hwTag, svtid_map_); 

        for (int ss = 0; ss < 6; ss++)
        {
            //histokey = hwTag + "_SvtHybrids_s"+std::to_string(ss)+"_hh";
            //histokey = hwTag + "_SvtHybrids_s0_hh";
            Fill2DHisto(histokey, 
                    (float)rawSvtHit->getStrip(),
                    (float)rawSvtHit->getADCs()[ss], 
                    weight);
        }
    }

    Event_number++;
}      
