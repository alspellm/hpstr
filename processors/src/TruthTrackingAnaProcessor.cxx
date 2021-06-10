#include "TruthTrackingAnaProcessor.h"
#include <iomanip>
#include "utilities.h"

TruthTrackingAnaProcessor::TruthTrackingAnaProcessor(const std::string& name, Process& process)
    : Processor(name, process) { 
    }

TruthTrackingAnaProcessor::~TruthTrackingAnaProcessor() { 
}

void TruthTrackingAnaProcessor::configure(const ParameterSet& parameters) {

    std::cout << "Configuring TruthTrackingAnaProcessor" << std::endl;
    try
    {
        debug_                = parameters.getInteger("debug",debug_);
        trkCollName_          = parameters.getString("trkCollName",trkCollName_);
        histCfgFilename_      = parameters.getString("histCfg",histCfgFilename_);
        doTruth_              = (bool) parameters.getInteger("doTruth",doTruth_);
        truthHistCfgFilename_ = parameters.getString("truthHistCfg",truthHistCfgFilename_);
        selectionCfg_         = parameters.getString("selectionjson",selectionCfg_); 
        purityCut_           = parameters.getDouble("puritycut", purityCut_);
        regionSelections_     = parameters.getVString("regionDefinitions",regionSelections_);
    }
    catch (std::runtime_error& error)
    {
        std::cout << error.what() << std::endl;
    }

}

void TruthTrackingAnaProcessor::initialize(TTree* tree) {

    //Init histos
    trkHistos_ = new TrackHistos(trkCollName_);
    trkHistos_->loadHistoConfig(histCfgFilename_);
    trkHistos_->doTrackComparisonPlots(false);
    trkHistos_->DefineHistos();
    // Init tree
    tree->SetBranchAddress(trkCollName_.c_str(), &tracks_, &btracks_);
    
    if (!selectionCfg_.empty()) {
        trkSelector_ = std::make_shared<BaseSelector>(name_+"_trkSelector",selectionCfg_);
        trkSelector_->setDebug(debug_);
        trkSelector_->LoadSelection();
    }
    
    if (doTruth_) {
        truthHistos_ = new TrackHistos(trkCollName_+"_truthComparison");
        truthHistos_->loadHistoConfig(histCfgFilename_);
        truthHistos_->DefineHistos();
        truthHistos_->loadHistoConfig(truthHistCfgFilename_);
        truthHistos_->DefineHistos();
        truthHistos_->doTrackComparisonPlots(false);
        
        //tree->SetBranchAddress(truthCollName_.c_str(),&truth_tracks_,&btruth_tracks_);
    }

    //Setup regions
    for (unsigned int i_reg = 0; i_reg < regionSelections_.size(); i_reg++){
        std::string regname = AnaHelpers::getFileName(regionSelections_[i_reg], false);
        std::cout << "Setting up region :: " << regname << std::endl;
        reg_selectors_[regname] = std::make_shared<BaseSelector>(regname, regionSelections_[i_reg]);
        reg_selectors_[regname]->setDebug(debug_);
        reg_selectors_[regname]->LoadSelection();
        
        reg_trackHistos_[regname] = std::make_shared<TrackHistos>(regname);
        reg_trackHistos_[regname]->loadHistoConfig(histCfgFilename_);
        reg_trackHistos_[regname]->doTrackComparisonPlots(false);
        reg_trackHistos_[regname]->DefineHistos();

        reg_truthHistos_[regname] = std::make_shared<TrackHistos>(regname +"_truthComparison");
        reg_truthHistos_[regname]->loadHistoConfig(histCfgFilename_);
        reg_truthHistos_[regname]->doTrackComparisonPlots(false);
        reg_truthHistos_[regname]->DefineHistos();

        regions_.push_back(regname);
    }

}

bool TruthTrackingAnaProcessor::process(IEvent* ievent) {

    double weight = 1.;
    // Loop over all the LCIO Tracks and add them to the HPS event.
    int n_sel_tracks = 0;
    for (int itrack = 0; itrack < tracks_->size(); ++itrack) {
        
        if (trkSelector_) trkSelector_->getCutFlowHisto()->Fill(0.,weight);
        
        // Get a track
        Track* track = tracks_->at(itrack);
        double purity = track->getTrackTruthPurity();
        if (purity < purityCut_)
            continue;
        int n2dhits_onTrack = !track->isKalmanTrack() ? track->getTrackerHitCount() * 2 : track->getTrackerHitCount();
        
        //Track Selection
        if (trkSelector_ && !trkSelector_->passCutGt("n_hits_gt",n2dhits_onTrack,weight))
            continue;
        
        if (trkSelector_ && !trkSelector_->passCutGt("pt_gt",fabs(track->getPt()),weight))
            continue;

        Track* truth_track = nullptr;

        //Get the truth track
        if (doTruth_) { 
            truth_track = (Track*) track->getTruthLink().GetObject();
            if (!truth_track)
                std::cout<<"Warnings::TruthTrackingAnaProcessor::Requested Truth track but couldn't find it in the ntuple"<<std::endl;
        }
        
        if(debug_ > 0)
        {
            std::cout<<"========================================="<<std::endl;
            std::cout<<"========================================="<<std::endl;
            std::cout<<"Track params:           "<<std::endl;
            track->Print();
        }
        
        trkHistos_->Fill1DHistograms(track);
        trkHistos_->Fill2DTrack(track);
        
        if (truthHistos_) {
            truthHistos_->Fill1DHistograms(truth_track);
            truthHistos_->Fill2DTrack(track);
            truthHistos_->Fill1DTrackTruth(track, truth_track);
        }

        //generate track hit code
        int hitCode = 0;
        int* goodhits = track->getTrackTruthGoodHits();
        for (int layer = 0; layer < 4; layer++) {
            if (goodhits[layer] != 0){
                std::cout << "hitCode: " << hitCode << std::endl;
                hitCode = hitCode | (0x1 << layer);
            }
        }

        std::string trkname = "";
        double charge = (double) track->getCharge();
        if (charge < 0)
            trkname = "ele_";
        else
            trkname = "pos_";

        for (auto region : regions_){
            if(debug_) std::cout << "Check for region " << region 
                << " hc " << hitCode
                << " lt:" << !reg_selectors_[region]->passCutLt("hitCode_lt", ((double)hitCode)-0.5, weight)
                << std::endl;
            //Hit code req
            if (!reg_selectors_[region]->passCutLt("hitCode_lt", ((double)hitCode)-0.5, weight) ) continue;
            if (!reg_selectors_[region]->passCutGt("hitCode_gt", ((double)hitCode)+0.5, weight) ) continue;
            reg_trackHistos_[region]->Fill1DTrack(track, weight, trkname);
            reg_truthHistos_[region]->Fill1DTrackTruth(track, truth_track, weight, trkname);
        }

        
        n_sel_tracks++;
    }//Loop on tracks

    trkHistos_->Fill1DHisto("n_tracks_h",n_sel_tracks);
    

    return true;
}

void TruthTrackingAnaProcessor::finalize() { 

    trkHistos_->saveHistos(outF_,trkCollName_);
    delete trkHistos_;
    trkHistos_ = nullptr;
    if (trkSelector_)
        trkSelector_->getCutFlowHisto()->Write();

    if (truthHistos_) {
        truthHistos_->saveHistos(outF_,trkCollName_+"_truth");
        delete truthHistos_;
        truthHistos_ = nullptr;
    }

    for (reg_it it = reg_trackHistos_.begin(); it!=reg_trackHistos_.end(); ++it){
        std::string dirName = it->first;
        (it->second)->saveHistos(outF_, dirName);
    }

    for (reg_it it = reg_truthHistos_.begin(); it!=reg_truthHistos_.end(); ++it){
        std::string dirName = it->first;
        (it->second)->saveHistos(outF_, dirName);
    }

    //trkHistos_->Clear();
}

DECLARE_PROCESSOR(TruthTrackingAnaProcessor); 
