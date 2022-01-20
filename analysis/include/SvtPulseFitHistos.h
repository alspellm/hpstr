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

        double GetHitTime(int sample_number, int cdel);
        std::map<std::string, TProfile*> getTProfiles(){ return tprofiles_;};
        double getAmplitudeIntegralNorm(double tau1, double tau2);
        TF1* fourPoleFitFunction();
        void saveHistos(TFile* outFile);
        void definePulseHistos(std::string name);
        void setOutFile(TFile* outFile){ outFile_ = outFile;};
        void initHistos();
        void jlab2019CalPulseScan(TTree* rawhitsTree);
        void adjustClock25nsTo24ns();
        void fitTGraphPulses();
        void passFitTupleOut(FlatTupleMaker* rawhitfits_tup){rawhitfits_tup_ = rawhitfits_tup;};
        void fit2DHistoPulses();
        void setSelectCalgroup(int calgroup){select_calgroup_ = calgroup;};
    

    private:

        int Event_number=0;
        int debug_ = 1;

        TH1F* svtCondHisto{nullptr};  
        
        //output fit tuple
        FlatTupleMaker* rawhitfits_tup_{nullptr};

        int select_calgroup_ = -1;

        //ModuleMapper
        ModuleMapper* mmapper_;
        std::map<std::string, std::map<int,int>> svtid_map_;
        TTree* hittree_{nullptr};
        RawSvtHit* rawhit_{nullptr};
        std::map<std::string,TProfile*> tprofiles_;
        std::map<std::string,TH2F*> histos2d_;
        std::map<std::string,TH1F*> histos1d_;

        //Pulse histograms
        std::map<std::string,TGraphErrors*> tgrapherrs_;
        std::map<std::string,TH2F*>pulsehistos2d_;
        std::map<std::string,TH1F*>baselines_;

        TFile* outFile_{nullptr};
        int nan_channels_ = 0;

    protected:
        typedef std::map<std::string, TProfile*>::iterator ittp;
        typedef std::map<std::string, TH1F*>::iterator it1d;
        typedef std::map<std::string, TH2F*>::iterator it2d;
        typedef std::map<std::string, TGraphErrors*>::iterator it_tgr;

};


#endif
