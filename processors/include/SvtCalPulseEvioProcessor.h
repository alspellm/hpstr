#ifndef __SVTBL2D_EVIOPROCESSOR_H__
#define __SVTBL2D_EVIOPROCESSOR_H__

//HPSTR
#include "Processor.h"
#include "ModuleMapper.h"
#include "RawSvtHit.h"
#include "SvtPulseFitHistos.h"
#include "FlatTupleMaker.h"
//EVIO
#include "EvioTool/HPSEvioReader.h"

//ROOT
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"

//CPLUSPLUS
#include <map>
#include <vector>
#include <memory>
#include <ratio>
#include <chrono>
#include <limits.h>
#include <iostream>
#include <string>

class SvtCalPulseEvioProcessor : public Processor {

    public:

        SvtCalPulseEvioProcessor(const std::string& name, Process& process);

        ~SvtCalPulseEvioProcessor();

        virtual void configure(const ParameterSet& parameters);

        virtual void initialize(std::string inFilename, std::string outFilename);

        virtual bool process();

        virtual void initialize(TTree* tree) {};

        virtual bool process(IEvent* event) { return true;};

        virtual void finalize();

    private:

        //Debug Level
        int debug_{0};
        SvtPulseFitHistos* svtPulseFitHistos{nullptr};
        int select_calgroup_;

        //Initialize some containers
        ModuleMapper * mmapper_;
        std::vector<RawSvtHit*> rawSvtHits_;
        TFile* outFile_{nullptr};
        
        //configuration parameters
        std::string histCfgFilename_;
        std::string histNames_{""};
        std::string chNumCfg_{""};
        std::string trigFilename_{""};
        std::string inFilename_{""};
        int processEvio_{1};
        int fitPulses_{1};
        int buildPulseHistos_{1};

        //evio interface 
        HPSEvioReader * etool{nullptr};
        
        int run_number_{0};
        unsigned long trigtime_start_{0};

        //Rawhit tuple
        FlatTupleMaker* rawhits_tup_{nullptr};
        //Pulse fit tuple
        FlatTupleMaker* rawhitfits_tup_{nullptr};

        //histos
        std::map<std::string, TH1F*> histos1d_;

        //tree
        TTree* rawhitsTree_{nullptr};
        //std::vector<RawSvtHit*> rawsvthits_;
        int eventnumber_;
};


#endif
