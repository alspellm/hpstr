/**
 * @file SvtCalPulseEvioProcessor.cxx
 * @author Emrys Peets, Stanford, SLAC National Accelerator Laboratory
 * @author Cameron Bravo, SLAC National Accelerator Laboratory
 */

#include "SvtCalPulseEvioProcessor.h"

SvtCalPulseEvioProcessor::SvtCalPulseEvioProcessor(const std::string& name, Process& process)
    : Processor(name, process) { 
        mmapper_ = new ModuleMapper();
    }

SvtCalPulseEvioProcessor::~SvtCalPulseEvioProcessor() { 
}

void SvtCalPulseEvioProcessor::configure(const ParameterSet& parameters) {

    std::cout << "Configuring SvtCalPulseEvioProcessor" << std::endl;
    try
    {
        debug_          = parameters.getInteger("debug");
        chNumCfg_   = parameters.getString("chNumCfg");
        trigFilename_   = parameters.getString("trigConf");
        histNames_  = parameters.getString("histNames");
        histCfgFilename_  = parameters.getString("histCfg");
    }
    catch (std::runtime_error& error)
    {
        std::cout << error.what() << std::endl;
    }
    if (histNames_ != "fw" && histNames_ != "sw") 
    {
        throw std::runtime_error("[ SvtCalPulseEvioProcessor ]: chNumCfg must be 'fw' or 'sw'!");
    }
    if (chNumCfg_ != "fw" && chNumCfg_ != "sw") 
    {
        throw std::runtime_error("[ SvtCalPulseEvioProcessor ]: chNumCfg must be 'fw' or 'sw'!");
    }
}

void SvtCalPulseEvioProcessor::initialize(std::string inFilename, std::string outFilename) {
    std::cout << "SvtCalPulseEvioProcessor::initialize" << std::endl;
    inFilename_ = inFilename;
    outF_ = new TFile(outFilename.c_str(),"RECREATE");
    etool = new HPSEvioReader();
    etool->Open(inFilename.c_str());
    etool->SVT->fSaveHeaders = true;

    if(trigFilename_.size()>2){
        etool->TrigConf->Parse_trigger_file(trigFilename_); // If requested, parse the Trigger config file supplied.
        etool->ECAL->Config();
        std::cout << "Parsed trigger file: " << trigFilename_ << std::endl;
    }

    //DEBUGGING  
    if(debug_==0){
        etool->fDebug = 0b000000;
    }else if(debug_==1){
        etool->fDebug = 0b000001;
    }else if(debug_==2){
        etool->fDebug = 0b000011;
    } else if(debug_ == 3){
        etool->fDebug = 0b000111;
    } else{
        etool->fDebug = 0xFF;
    }


    bool found = false;
    while( !found &&   etool->Next()== S_SUCCESS){
        if( (etool->this_tag & 128) == 128){
            found = true;
            run_number_ = etool->GetRunNumber();
            trigtime_start_ = etool->GetTrigTime();
            break;
        }
    }
    if( found == false) std::cout << "WARNING -- Not able to find a bank with a runnumber! \n";
    etool->Close();

    if(etool->SVT){
        if(run_number_ < 8000){ // This is 2015 or 2016 data.
            etool->Set2016Data();
        }else{
            etool->Set2019Data();
        }
    }else{
        std::cout << "NO SVT initialized \n";
    }

    etool->fAutoAdd = false;

    if(debug_>1) etool->PrintBank(5);

    //Setup flat tuple branches
    //rawhits_tup_ = new FlatTupleMaker(outFilename.c_str(), "rawhits");
    rawhits_tup_ = new FlatTupleMaker("rawhits");
    rawhits_tup_->setOutFile(outF_);
    //hardware tag F<n>H<m>
    rawhits_tup_->addString("hwTag");
    rawhits_tup_->addVariable("layer");
    rawhits_tup_->addVariable("module");
    rawhits_tup_->addVariable("channel");
    rawhits_tup_->addVariable("svtid");
    rawhits_tup_->addVariable("cdel");
    rawhits_tup_->addVariable("calgroup");
    rawhits_tup_->addVariable("adc0");
    rawhits_tup_->addVariable("adc1");
    rawhits_tup_->addVariable("adc2");
    rawhits_tup_->addVariable("adc3");
    rawhits_tup_->addVariable("adc4");
    rawhits_tup_->addVariable("adc5");
    rawhits_tup_->addVariable("event");

    //tuple for storing rawhit fits
    //rawhitfits_tup_ = new FlatTupleMaker(outFilename.c_str(), "fits");
    rawhitfits_tup_ = new FlatTupleMaker("fits");
    rawhitfits_tup_->setOutFile(outF_);
    rawhitfits_tup_->addVariable("svtid");
    rawhitfits_tup_->addVariable("layer");
    rawhitfits_tup_->addVariable("module");
    rawhitfits_tup_->addVariable("t0");
    rawhitfits_tup_->addVariable("tau1");
    rawhitfits_tup_->addVariable("tau2");
    rawhitfits_tup_->addVariable("amp");
    rawhitfits_tup_->addVariable("t0err");
    rawhitfits_tup_->addVariable("tau1err");
    rawhitfits_tup_->addVariable("tau2err");
    rawhitfits_tup_->addVariable("amperr");
    rawhitfits_tup_->addVariable("integralNorm");
    rawhitfits_tup_->addVariable("baseline");
    rawhitfits_tup_->addVariable("chi2");
    rawhitfits_tup_->addVariable("ndf");

    //Init histos
    svtPulseFitHistos = new SvtPulseFitHistos("raw_hits", mmapper_);
}

bool SvtCalPulseEvioProcessor::process() {

    int maxevents = 50;
    int eventn = 0;

    std::cout << "SvtCalPulseEvioProcessor::process" << std::endl;
    unsigned long evt_count=0;
    int l0APVmap[4] = {1, 0, 2, 3};
    int APVmap[5] = {4, 3, 2, 1, 0};

    etool->Open(inFilename_.c_str());
    std::cout << "SvtCalPulseEvioProcessor::opened file" << std::endl;
    while(etool->Next() == S_SUCCESS){
        if(etool->Head->GetEventNumber() == 1 || (etool->Head->GetEventNumber())%2==0)
            continue;
        if( (etool->this_tag & 128) != 128) continue;
        if(debug_) cout<<"EVIO Event " << etool->Head->GetEventNumber() << endl;
        if(debug_) cout << "Event Number:  " << etool->Head->GetEventNumber() << "  seq: " << evt_count << endl;
        cout<<"EVIO Event " << etool->Head->GetEventNumber() << endl;
        eventn++;
        if(eventn > maxevents){
            break;
        }
        for(int i = 0; i < rawSvtHits_.size(); i++)
        {
            delete rawSvtHits_[i];
        }
        rawSvtHits_.clear();
        //      etool->VtpTop->ParseBank();
        //      etool->VtpBot->ParseBank();
        rawhits_tup_->setVariableValue("event", (double)etool->Head->GetEventNumber());
        evt_count++;

        if(debug_>0) {
            etool->PrintBank(10);
        }
        if(etool->SVT){
            //etool->PrintBank(0);
            //std::cout << etool->SVT->svt_data.size() << " raw svt hits" << endl;
            for(int i = 0; i < etool->SVT->svt_data.size(); i++)
            {
                RawSvtHit * rawHit = new RawSvtHit();
                std::string hwTag = "F"+std::to_string(int(etool->SVT->svt_data[i].head.feb_id))+"H"+std::to_string(int(etool->SVT->svt_data[i].head.hyb_id));
                std::string swTag = mmapper_->getSwFromHw(hwTag);
                int layer = std::stoi(swTag.substr(2, swTag.find("_")));
                int module = std::stoi(swTag.substr(swTag.find("m")+1, swTag.length()));
                rawHit->setLayer(layer);
                rawHit->setModule(module);
                int apv = int(etool->SVT->svt_data[i].head.apv);
                int strip = apv*128;
                if (chNumCfg_ == "sw")
                {
                    strip = l0APVmap[apv]*128;
                    if ( int(etool->SVT->svt_data[i].head.feb_id) > 1 ) strip = APVmap[apv]*128;
                }
                strip += int(etool->SVT->svt_data[i].head.chan);
                rawHit->setStrip(strip);
                int adcs[6];
                for(int adcI = 0; adcI < etool->SVT->svt_data.size(); adcI++)
                {
                    adcs[adcI] = int(etool->SVT->svt_data[i].samples[adcI]);
                }
                rawHit->setADCs(adcs);
                rawSvtHits_.push_back(rawHit);
            }

            svtPulseFitHistos->buildRawSvtHitsTuple(&rawSvtHits_,rawhits_tup_);
        }

    }

    return true;
}

void SvtCalPulseEvioProcessor::finalize() { 

    std::cout << "SvtCalPulseEvioProcessor::finalize" << std::endl;
    outF_->cd();
    std::cout << "SvtCalPulseEvioProcessor::write rawhits tuple" << std::endl;
    /*
    rawhits_tup_->close();
    outF_->Close();
    */
    rawhits_tup_->writeTree(outF_);
    
    /*
    //Read Ttree and pass to fit pulses
    outF_ = new TFile("testout.root","UPDATE");
    outF_->ls();
    */
    TTree* rawhittree = (TTree*)outF_->Get("rawhits");
    rawhittree->Print();
    std::cout << "Fitting raw hit pulses" << std::endl;
    svtPulseFitHistos->fitRawHitPulses(rawhittree, rawhitfits_tup_);
    std::cout << "save fit tuple" << std::endl;
    //outF_->cd();
    //rawhitfits_tup_->writeTree();
    rawhitfits_tup_->writeTree(outF_);
    svtPulseFitHistos->saveTProfiles(outF_);
    svtPulseFitHistos->saveHistos(outF_);
    
    //svtPulseFitHistos->saveHistos(outF_,"");
    outF_->Close();
    delete svtPulseFitHistos;
    svtPulseFitHistos = nullptr;
    delete mmapper_;
    mmapper_ = nullptr;

}

DECLARE_PROCESSOR(SvtCalPulseEvioProcessor); 
