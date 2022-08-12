#include "DebugProcessor.h"
#include <iomanip>
#include "utilities.h"
#include "Particle.h"
#include "CalCluster.h"

DebugProcessor::DebugProcessor(const std::string& name, Process& process)
    : Processor(name, process) { 
    }

DebugProcessor::~DebugProcessor() { 
}

void DebugProcessor::configure(const ParameterSet& parameters) {

    std::cout << "Configuring DebugProcessor" << std::endl;
    try
    {
        debug_                = parameters.getInteger("debug",debug_);
        trkCollName_          = parameters.getString("trkCollName",trkCollName_);
        vtxCollName_          = parameters.getString("vtxCollName",vtxCollName_);
        histCfgFilename_      = parameters.getString("histCfg",histCfgFilename_);
        doTruth_              = (bool) parameters.getInteger("doTruth",doTruth_);
        truthHistCfgFilename_ = parameters.getString("truthHistCfg",truthHistCfgFilename_);
        selectionCfg_         = parameters.getString("selectionjson",selectionCfg_); 
        isData_               = parameters.getInteger("isData",isData_);

        //region definitions
        regionSelections_ = parameters.getVString("regionDefinitions",regionSelections_);
    }
    catch (std::runtime_error& error)
    {
        std::cout << error.what() << std::endl;
    }

}

void DebugProcessor::initialize(TTree* tree) {

    //Init histos
    trkHistos_ = new TrackHistos(trkCollName_);
    trkHistos_->loadHistoConfig(histCfgFilename_);
    trkHistos_->doTrackComparisonPlots(false);
    trkHistos_->DefineHistos();
    trkHistos_->DefineTrkHitHistos();
    // Init tree
    tree->SetBranchAddress("EventHeader", &evth_ , &bevth_);
    tree->SetBranchAddress(trkCollName_.c_str(), &tracks_, &btracks_);
    tree->SetBranchAddress(vtxCollName_.c_str(), &vtxs_ , &bvtxs_);
    
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

}

bool DebugProcessor::process(IEvent* ievent) {

    double weight = 1.;
    // Loop over all the LCIO Tracks and add them to the HPS event.
    int n_sel_tracks = 0;
    for (int itrack = 0; itrack < tracks_->size(); ++itrack) {
        
        if (trkSelector_) trkSelector_->getCutFlowHisto()->Fill(0.,weight);
        
        // Get a track
        Track* track = tracks_->at(itrack);
        int charge = track->getCharge();
        int n2dhits_onTrack = !track->isKalmanTrack() ? track->getTrackerHitCount() * 2 : track->getTrackerHitCount();
        double momentum = track->getP();
        double trackTime = track->getTrackTime();
        double tanlambda = track->getTanLambda();
        bool isTop = track->isTopTrack();
        bool isPos = false;
        if (charge > 0 )
            isPos = true;
        /*
        for ( int i_vtx = 0; i_vtx <  vtxs_->size(); i_vtx++ ) {
            Vertex* vtx = vtxs_->at(i_vtx);
            Particle* ele = nullptr;
            Track* ele_trk = nullptr;
            Particle* pos = nullptr;
            Track* pos_trk = nullptr;
        }*/

        //Debug tracking plots

        //Track Selection

        if (trkSelector_ && !trkSelector_->passCutGt("n_hits_gt",n2dhits_onTrack,weight))
            continue;

        if (trkSelector_ && !trkSelector_->passCutLt("p_lt",fabs(track->getP()),weight))
            continue;
        
        if (trkSelector_ && !trkSelector_->passCutGt("p_gt",fabs(track->getPt()),weight))
            continue;

        //Check data event for Trigger
        if (isData_) {
            if (trkSelector_ && !trkSelector_->passCutEq("Pair1_eq",(int)evth_->isPair1Trigger(),weight))
                continue;
        }
        //Ele Track Quality - Chi2Ndf
        if (!trkSelector_->passCutLt("chi2ndf_lt",track->getChi2Ndf(),weight))
            continue;
        

        Track* truth_track = nullptr;

        //Get the truth track
        if (doTruth_) { 
            truth_track = (Track*) track->getTruthLink().GetObject();
            if (!truth_track)
                std::cout<<"Warnings::DebugProcessor::Requested Truth track but couldn't find it in the ntuple"<<std::endl;
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

        if(isTop&&isPos) trkHistos_->Fill1DTrack(track, weight, "topPos_");
        if(isTop&&!isPos) trkHistos_->Fill1DTrack(track, weight, "topEle_");
        if(!isTop&&isPos) trkHistos_->Fill1DTrack(track, weight, "botPos_");
        if(!isTop&&!isPos) trkHistos_->Fill1DTrack(track, weight, "botEle_");

        if(isTop&&isPos) trkHistos_->Fill2DTrack(track, weight, "topPos_");
        if(isTop&&!isPos) trkHistos_->Fill2DTrack(track, weight, "topEle_");
        if(!isTop&&isPos) trkHistos_->Fill2DTrack(track, weight, "botPos_");
        if(!isTop&&!isPos) trkHistos_->Fill2DTrack(track, weight, "botEle_");

        //nhits plots
        if(isTop&&isPos) trkHistos_->Fill1DTrack(track, weight, "topPos_nhits_"+std::to_string(n2dhits_onTrack)+"_");
        if(isTop&&!isPos) trkHistos_->Fill1DTrack(track, weight, "topEle_nhits_"+std::to_string(n2dhits_onTrack)+"_");
        if(!isTop&&isPos) trkHistos_->Fill1DTrack(track, weight, "botPos_nhits_"+std::to_string(n2dhits_onTrack)+"_");
        if(!isTop&&!isPos) trkHistos_->Fill1DTrack(track, weight, "botEle_nhits_"+std::to_string(n2dhits_onTrack)+"_");

        if(isTop&&isPos) trkHistos_->Fill2DTrack(track, weight, "topPos_nhits_"+std::to_string(n2dhits_onTrack)+"_");
        if(isTop&&!isPos) trkHistos_->Fill2DTrack(track, weight, "topEle_nhits_"+std::to_string(n2dhits_onTrack)+"_");
        if(!isTop&&isPos) trkHistos_->Fill2DTrack(track, weight, "botPos_nhits_"+std::to_string(n2dhits_onTrack)+"_");
        if(!isTop&&!isPos) trkHistos_->Fill2DTrack(track, weight, "botEle_nhits_"+std::to_string(n2dhits_onTrack)+"_");
        
        if (truthHistos_) {
            truthHistos_->Fill1DHistograms(truth_track);
            truthHistos_->Fill2DTrack(track);
        }
        
        n_sel_tracks++;
    }//Loop on tracks

    trkHistos_->Fill1DHisto("n_tracks_h",n_sel_tracks);

    return true;
}

void DebugProcessor::finalize() { 

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
    //trkHistos_->Clear();
}

DECLARE_PROCESSOR(DebugProcessor); 
