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
        //histCfgFilename_  = parameters.getString("histCfg");
        processEvio_      = parameters.getInteger("processEvio");
        fitPulses_       = parameters.getInteger("fitPulses");
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

    std::cout << "processEvio? " << processEvio_ << std::endl;
    std::cout << "fitPulses? " << fitPulses_ << std::endl;
}

void SvtCalPulseEvioProcessor::initialize(std::string inFilename, std::string outFilename) {

    std::cout << "SvtCalPulseEvioProcessor::initialize" << std::endl;
    inFilename_ = inFilename;
    outFile_ = new TFile(outFilename.c_str(),"RECREATE");
    if(processEvio_){
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
    }
    rawhitsTree_ = new TTree("rawsvthits","rawsvthits");
    rawhitsTree_->Branch("event",&eventnumber_);
    rawhitsTree_->Branch("rawsvthits",&rawSvtHits_);

    //Define tuple to store channel pulse fit results
    rawhitfits_tup_ = new FlatTupleMaker("fits");
    rawhitfits_tup_->setOutFile(outFile_);
    rawhitfits_tup_->addString("hwTag");
    rawhitfits_tup_->addVariable("svtid");
    rawhitfits_tup_->addVariable("channel");
    rawhitfits_tup_->addVariable("layer");
    rawhitfits_tup_->addVariable("module");
    rawhitfits_tup_->addVariable("t0");
    rawhitfits_tup_->addVariable("tau1");
    rawhitfits_tup_->addVariable("tau2");
    rawhitfits_tup_->addVariable("amp");
    rawhitfits_tup_->addVariable("baseline");
    rawhitfits_tup_->addVariable("chi2");
    rawhitfits_tup_->addVariable("ndf");
    rawhitfits_tup_->addVariable("t0err");
    rawhitfits_tup_->addVariable("tau1err");
    rawhitfits_tup_->addVariable("tau2err");
    rawhitfits_tup_->addVariable("amperr");
    rawhitfits_tup_->addVariable("integralNorm");
    rawhitfits_tup_->addVariable("nanfit");

    //Init histos
    if(fitPulses_){
        svtPulseFitHistos = new SvtPulseFitHistos("raw_hits", mmapper_);
        svtPulseFitHistos->passFitTupleOut(rawhitfits_tup_);
        svtPulseFitHistos->initHistos();
    }
}

bool SvtCalPulseEvioProcessor::process() {

    if(!processEvio_){
        std::cout << "SKIPPING EVIO PROCESSOR" << std::endl;
        return true;
    }
    //Else input is evio data, need to process evio data into a TTree before fitting
    else {
        std::cout << "PROCESSING EVIO DATA INTO TTREE" << std::endl;
        int eventn = 0;
        int maxevents = 20;
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
            eventn++;
            //if (eventn > maxevents){
            //    break;
            //}
            for(int i = 0; i < rawSvtHits_.size(); i++)
            {
                delete rawSvtHits_[i];
            }
            rawSvtHits_.clear();
            //      etool->VtpTop->ParseBank();
            //      etool->VtpBot->ParseBank();
            //rawsvthits_tup_->setVariableValue("event", (double)etool->Head->GetEventNumber());
            eventnumber_ = (double)etool->Head->GetEventNumber();
            //std::cout << eventnumber_ << std::endl;

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
                    //rawsvthits_.push_back(rawHit);
                }

                
                //svtPulseFitHistos->buildRawSvtHitsTuple(&rawSvtHits_,rawsvthits_tup_);
            }
            outFile_->cd(); 
            rawhitsTree_->Fill();
        }
        return true;
    }
}

void SvtCalPulseEvioProcessor::finalize() { 

    std::cout << "SvtCalPulseEvioProcessor::finalize" << std::endl;

    TTree* rawhitTree{nullptr};
    //If input file is evio, write evio->rawhit TTree
    if (processEvio_){
        outFile_->cd();
        std::cout << "SAVING EVIO DATA TO TTREE" << std::endl;
        outFile_->Write();
    }

    if (processEvio_ && fitPulses_){
        outFile_->cd();
        std::cout << "READING PROCESSED EVIO TTREE TO FIT PULSES" << std::endl;
        rawhitTree = (TTree*)outFile_->Get("rawsvthits");
    }

    else if (fitPulses_ && !processEvio_){
        TFile* readF = new TFile(inFilename_.c_str(),"READ");
        readF->cd();
        std::cout << "READING TTREE FROM INPUT FILE TO FIT PULSES" << std::endl;
        rawhitTree = (TTree*)readF->Get("rawsvthits");
        std::cout << "Tree Read" << std::endl;
    }

    if(fitPulses_){
        std::cout << "FITTING CALIBRATION SCAN PULSES" << std::endl;
        svtPulseFitHistos->jlab2019CalPulseScan(rawhitTree);
        svtPulseFitHistos->saveHistos(outFile_);
    }

    outFile_->Close();
    delete svtPulseFitHistos;
    svtPulseFitHistos = nullptr;
    delete mmapper_;
    mmapper_ = nullptr;

}

DECLARE_PROCESSOR(SvtCalPulseEvioProcessor); 
