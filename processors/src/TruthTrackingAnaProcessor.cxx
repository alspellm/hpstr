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
    
    truthHistos_ = new TrackHistos(trkCollName_+"_truthComparison");
    truthHistos_->loadHistoConfig(histCfgFilename_);
    truthHistos_->DefineHistos();
    truthHistos_->loadHistoConfig(truthHistCfgFilename_);
    truthHistos_->DefineHistos();
    truthHistos_->doTrackComparisonPlots(false);
        
        //tree->SetBranchAddress(truthCollName_.c_str(),&truth_tracks_,&btruth_tracks_);

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
        reg_truthHistos_[regname]->loadHistoConfig(truthHistCfgFilename_);
        reg_truthHistos_[regname]->doTrackComparisonPlots(false);
        reg_truthHistos_[regname]->DefineHistos();

        reg_misHistos_[regname] = std::make_shared<TrackHistos>(regname +"_misassHits");
        reg_misHistos_[regname]->loadHistoConfig(truthHistCfgFilename_);
        reg_misHistos_[regname]->doTrackComparisonPlots(false);
        reg_misHistos_[regname]->DefineHistos();

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
        std::string trkname = "";
        bool isPos = (track->getOmega() < 0.0);
        int trkType = (int)isPos;
        double charge = (double) track->getCharge();
        if (charge < 0)
            trkname = "ele_";
        else
            trkname = "pos_";
        double purity = track->getTrackTruthPurity();
        if (purity < purityCut_)
            continue;
        //adding cuts for debug
        if (track->getP() < 4)
            continue;
        int n2dhits_onTrack = !track->isKalmanTrack() ? track->getTrackerHitCount() * 2 : track->getTrackerHitCount();
        //added cut for debug
        if (n2dhits_onTrack < 14)
            continue;
        
        //Track Selection
        if (trkSelector_ && !trkSelector_->passCutGt("n_hits_gt",n2dhits_onTrack,weight))
            continue;
        
        if (trkSelector_ && !trkSelector_->passCutGt("pt_gt",fabs(track->getPt()),weight))
            continue;

        Track* truth_track = nullptr;

        //trkHistos_->Fill1DHistograms(track);
        trkHistos_->Fill1DTrack(track, weight, trkname);
        trkHistos_->Fill2DTrack(track);
        
        //Get the truth track
        truth_track = (Track*) track->getTruthLink().GetObject();
        //if no truth track, continue to next track
        if(!truth_track)
            continue;

        //truthHistos_->Fill1DHistograms(truth_track, weight, trkname);
        truthHistos_->Fill1DTrack(truth_track, weight, trkname);
        truthHistos_->Fill2DTrack(truth_track,weight, trkname);
        truthHistos_->Fill1DTrackTruth(track, truth_track, weight, trkname);

        //generate track hit code
        int hitCode = 0;
        int mishitCode = 0;
        int* goodhits = track->getTrackTruthGoodHits();
        for (int layer = 0; layer < 4; layer++) {
            if (goodhits[layer] != 0){
                hitCode = hitCode | (0x1 << layer);
            }
        }
        //Check hit misassociation for Tracks that have hits in first 4 layers
        if (hitCode == 15){
            for (int layer = 0; layer < 4; layer++) {
                if(goodhits[layer] == -1){
                    mishitCode = mishitCode | (0x1 << layer);
                }
            }
        }

        truthHistos_->Fill1DHisto("hitCode_h", hitCode);
        truthHistos_->Fill2DHisto("hitCode_trkType_hh", hitCode, trkType);
        truthHistos_->Fill1DHisto("mishitCode_h", mishitCode);
        truthHistos_->Fill2DHisto("mishitCode_trkType_hh", mishitCode, trkType);

        for (auto region : regions_){
            //Hit code req
            if (!reg_selectors_[region]->passCutLt("hitCode_lt", ((double)hitCode)-0.5, weight) ) continue;
            if (!reg_selectors_[region]->passCutGt("hitCode_gt", ((double)hitCode)+0.5, weight) ) continue;
            reg_trackHistos_[region]->Fill1DTrack(track, weight, trkname);
            reg_truthHistos_[region]->Fill1DTrackTruth(track, truth_track, weight, trkname);

            if (!reg_selectors_[region]->passCutLt("hitCode_lt", ((double)mishitCode)-0.5, weight) ) continue;
            if (!reg_selectors_[region]->passCutGt("hitCode_gt", ((double)mishitCode)+0.5, weight) ) continue;
            reg_misHistos_[region]->Fill1DTrackTruth(track, truth_track, weight, trkname);
        }
        
        n_sel_tracks++;
    }//Loop on tracks

    //trkHistos_->Fill1DHisto("n_tracks_h",n_sel_tracks);
    

    return true;
}

void TruthTrackingAnaProcessor::finalize() { 

    trkHistos_->saveHistos(outF_,trkCollName_);
    delete trkHistos_;
    trkHistos_ = nullptr;
    if (trkSelector_)
        trkSelector_->getCutFlowHisto()->Write();

    truthHistos_->saveHistos(outF_,trkCollName_+"_truth");
    delete truthHistos_;
    truthHistos_ = nullptr;

    for (reg_it it = reg_trackHistos_.begin(); it!=reg_trackHistos_.end(); ++it){
        std::string dirName = it->first + "_trackHistos";
        (it->second)->saveHistos(outF_, dirName);
    }

    for (reg_it it = reg_truthHistos_.begin(); it!=reg_truthHistos_.end(); ++it){
        std::string dirName = it->first + "_truthHistos";
        (it->second)->saveHistos(outF_, dirName);
    }

    for (reg_it it = reg_misHistos_.begin(); it!=reg_misHistos_.end(); ++it){
        std::string dirName = it->first + "_mishitHistos";
        (it->second)->saveHistos(outF_, dirName);
    }

    //trkHistos_->Clear();
}

DECLARE_PROCESSOR(TruthTrackingAnaProcessor); 
