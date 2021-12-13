#ifndef SVTPULSEFITHISTOS_H
#define SVTPULSEFITHISTOS_H

#include "TFile.h"
#include "HistoManager.h"
#include "TGraphErrors.h"
#include "TKey.h"
#include "TList.h"
#include "TH1.h"
#include "TrackerHit.h"
#include "RawSvtHit.h"
#include "TTree.h"
#include "FlatTupleMaker.h"
#include "TProfile.h"
#include "TF1.h"

#include "ModuleMapper.h"

#include <string>


class SvtPulseFitHistos : public HistoManager{

    public:
        SvtPulseFitHistos(const std::string& inputName, ModuleMapper* mmapper_);
        ~SvtPulseFitHistos();

        void FillHistogramsByHw(std::vector<RawSvtHit*> *rawSvtHits_,float weight = 1.);
        double GetHitTime(int sample_number, int cdel);
        void buildRawSvtHitsTuple(std::vector<RawSvtHit*> *rawSvtHits_, FlatTupleMaker* rawhits_tup_);
        void fitRawHitPulses(TTree* rawhittree);
        void defineTProfiles(int maxchannels);
        void defineTProfile(std::string name);
        std::map<std::string, TProfile*> getTProfiles(){ return tprofiles_;};
        virtual void saveTProfiles(TFile* outF = nullptr,std::string folder = "");
        double getAmplitudeIntegralNorm(double tau1, double tau2);
        TF1* fourPoleFitFunction();
        void fitPulse(TProfile* tprofile);


    private:

        int Event_number=0;
        int debug_ = 1;

        TH1F* svtCondHisto{nullptr};  
        //ModuleMapper
        ModuleMapper* mmapper_;
        std::map<std::string, std::map<int,int>> svtid_map_;
        TTree* hittree_{nullptr};
        RawSvtHit* rawhit_{nullptr};
        std::map<std::string,TProfile*> tprofiles_;

    protected:
        typedef std::map<std::string, TProfile*>::iterator ittp;

};


#endif
