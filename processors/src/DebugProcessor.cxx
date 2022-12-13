#include "DebugProcessor.h"
#include <iomanip>
#include "utilities.h"

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
        ecalColl_             = parameters.getString("ecalColl",ecalColl_);
        histCfgFilename_      = parameters.getString("histCfg",histCfgFilename_);
        doTruth_              = (bool) parameters.getInteger("doTruth",doTruth_);
        truthHistCfgFilename_ = parameters.getString("truthHistCfg",truthHistCfgFilename_);
        selectionCfg_         = parameters.getString("selectionjson",selectionCfg_); 
        timeOffset_           = parameters.getDouble("CalTimeOffset",timeOffset_);
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
    tree->SetBranchAddress(ecalColl_.c_str(), &ecal_  , &becal_);
    
    _ah = std::make_shared<AnaHelpers>();

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
    trkHistos_->Fill1DHisto("events_h",1);

    for ( int i_ecal = 0; i_ecal < ecal_->size(); i_ecal++ ) {

        //Check data event for Trigger
        if (isData_) {
            if (trkSelector_ && !trkSelector_->passCutEq("Pair1_eq",(int)evth_->isPair1Trigger(),0.0))
                continue;
        }

        CalCluster* cluster = ecal_->at(i_ecal);
        double clusterEnergy = cluster->getEnergy();
        double clusterTime = cluster->getTime() - timeOffset_;

        trkHistos_->Fill1DHisto("cluster_energy_h",clusterEnergy);
        trkHistos_->Fill1DHisto("cluster_time_h",clusterTime);
    }
    
    /*
    //Loop over vertices
    for(int i_vtx = 0; i_vtx < vtxs_->size(); i_vtx++) {

        if (isData_) {
            if (trkSelector_ && !trkSelector_->passCutEq("Pair1_eq",(int)evth_->isPair1Trigger(),weight))
                continue;
        }

        Vertex* vtx = vtxs_->at(i_vtx);

        Particle* ele = nullptr;
        Particle* pos = nullptr;

        _ah->GetParticlesFromVtx(vtx,ele,pos);
        CalCluster eleClus = ele->getCluster();
        CalCluster posClus = pos->getCluster();
        double eleClusterEnergy = eleClus.getEnergy();
        double posClusterEnergy = posClus.getEnergy();

        //Compute analysis variables here.
        //
        double ele_E = ele->getEnergy();
        double pos_E = pos->getEnergy();

        Track ele_trk = ele->getTrack();
        Track pos_trk = pos->getTrack();

        //Get the shared info - TODO change and improve

        Track* ele_trk_gbl = nullptr;
        Track* pos_trk_gbl = nullptr;

        if (!trkCollName_.empty()) {
            bool foundTracks = _ah->MatchToGBLTracks(ele_trk.getID(),pos_trk.getID(),
                    ele_trk_gbl, pos_trk_gbl, *tracks_);

            if (!foundTracks) {
                if (debug_)
                    std::cout<<"VertexAnaProcessor::ERROR couldn't find ele/pos in the "<<trkCollName_ <<"collection"<<std::endl;
                continue;
            }
        }
        else {

            ele_trk_gbl = (Track*) ele_trk.Clone();
            pos_trk_gbl = (Track*) pos_trk.Clone();
        }
        TVector3 recEleP(ele->getMomentum()[0],ele->getMomentum()[1],ele->getMomentum()[2]);
        TLorentzVector p_ele;
        p_ele.SetPxPyPzE(ele_trk_gbl->getMomentum()[0],ele_trk_gbl->getMomentum()[1],ele_trk_gbl->getMomentum()[2], ele_E);
        TLorentzVector p_pos;
        p_pos.SetPxPyPzE(pos_trk_gbl->getMomentum()[0],pos_trk_gbl->getMomentum()[1],pos_trk_gbl->getMomentum()[2], pos_E);

        double ele_trk_time = ele_trk_gbl->getTrackTime();
        double pos_trk_time = pos_trk_gbl->getTrackTime();
        trkHistos_->Fill1DHisto("ele_trk_time_h", ele_trk_time, weight);
        trkHistos_->Fill1DHisto("pos_trk_time_h", pos_trk_time, weight);


        double corr_eleClusterTime = ele->getCluster().getTime() - timeOffset_;
        double corr_posClusterTime = pos->getCluster().getTime() - timeOffset_;

        double botClusTime = 0.0;
        if(ele->getCluster().getPosition().at(1) < 0.0) botClusTime = ele->getCluster().getTime();
        else botClusTime = pos->getCluster().getTime();

        trkHistos_->Fill1DHisto("corr_eleClusterTime_h",corr_eleClusterTime, weight);
        trkHistos_->Fill1DHisto("corr_posClusterTime_h",corr_posClusterTime,weight);

        trkHistos_->Fill1DHisto("ele_trk_clust_dt_h", ele_trk_time-corr_eleClusterTime, weight);
        trkHistos_->Fill1DHisto("pos_trk_clust_dt_h", pos_trk_time-corr_posClusterTime, weight);

        std::cout << "trk time: " << ele_trk_time << std::endl;
        std::cout << "cluster time: " << ele->getCluster().getTime() << std::endl;
        std::cout << "time offset: " << timeOffset_ << std::endl;
        std::cout << "corrected clus time: " << corr_eleClusterTime << std::endl;
        std::cout << "track clus time diff: " << ele_trk_time-corr_eleClusterTime << std::endl;


        //Loop over track hits
        for(int ihit =0; ihit < ele_trk_gbl->getSvtHits()->GetEntries(); ihit++){
            TrackerHit* hit = (TrackerHit*) ele_trk_gbl->getSvtHits()->At(ihit);
            double hittime = hit->getTime();
            int layer = hit->getLayer();
            trkHistos_->Fill1DHisto("ele_trk_hit_time_h", hittime, weight);
            trkHistos_->Fill2DHisto("ele_trk_hit_time_v_layer_hh", layer, hittime, weight);

            std::cout << "LENGTH " << hit->getRawHits()->GetEntries() << std::endl;
            for(int iraw=0; iraw < hit->getRawHits()->GetEntries(); iraw++){
                RawSvtHit* rawhit = (RawSvtHit*) hit->getRawHits()->At(iraw);
                std::cout << "have rawhit" << std::endl;
                if(rawhit != NULL){
                    int sensor = rawhit->getSensor();
                    std::cout << "SENSOR: " << sensor << std::endl;
                }
            }

        }
        for(int ihit =0; ihit < pos_trk_gbl->getSvtHits()->GetEntries(); ihit++){
            TrackerHit* hit = (TrackerHit*) pos_trk_gbl->getSvtHits()->At(ihit);
            double hittime = hit->getTime();
            int layer = hit->getLayer();
            trkHistos_->Fill1DHisto("pos_trk_hit_time_h", hittime, weight);
            trkHistos_->Fill2DHisto("pos_trk_hit_time_v_layer_hh", layer, hittime, weight);
        }


        double weight = 0.0;
        //trkHistos_->Fill1DTrack(ele_trk_gbl,weight,"ele_");       
        //trkHistos_->Fill1DTrack(pos_trk_gbl,weight,"pos_");       
        //trkHistos_->Fill2DTrack(ele_trk_gbl,weight,"ele_");       
        //trkHistos_->Fill2DTrack(pos_trk_gbl,weight,"pos_");       

    }*/

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
        std::string charge_str = "ele";
        if (isPos)
            charge_str = "pos";

        //Check data event for Trigger
        if (isData_) {
            if (trkSelector_ && !trkSelector_->passCutEq("Pair1_eq",(int)evth_->isPair1Trigger(),weight))
                continue;
        }

        if (trkSelector_ && !trkSelector_->passCutGt("n_hits_gt",n2dhits_onTrack,weight))
            continue;
         if (trkSelector_ && !trkSelector_->passCutLt("trkTime_lt",fabs(track->getTrackTime()),weight))
            continue;
        if (trkSelector_ && !trkSelector_->passCutLt("chi2ndf_lt",track->getChi2Ndf(),weight))
            continue;
        if (trkSelector_ && !trkSelector_->passCutLt("p_lt",fabs(track->getP()),weight))
            continue;
        if (trkSelector_ && !trkSelector_->passCutGt("p_gt",fabs(track->getPt()),weight))
            continue;

        double trk_time = track->getTrackTime();
        trkHistos_->Fill1DHisto(charge_str+"_trk_time_h", trk_time, weight);
        trkHistos_->Fill1DHisto(charge_str+"_trk_p_h", momentum, weight);
        trkHistos_->Fill2DHisto(charge_str+"_trk_time_v_p_hh", momentum, trk_time, weight);

        //Loop over track hits
        for(int ihit =0; ihit < track->getSvtHits()->GetEntries(); ihit++){
            TrackerHit* hit = (TrackerHit*) track->getSvtHits()->At(ihit);
            double hittime = hit->getTime();
            int layer = hit->getLayer();
            trkHistos_->Fill1DHisto(charge_str+"_trk_strip_cluster_time_h", hittime, weight);
            trkHistos_->Fill2DHisto(charge_str+"_trk_strip_cluster_time_v_layer_hh", layer, hittime, weight);

            for(int iraw=0; iraw < hit->getRawHits()->GetEntries(); iraw++){
                RawSvtHit* rawhit = (RawSvtHit*) hit->getRawHits()->At(iraw);
                if(rawhit != NULL){
                    int sensor = rawhit->getSensor();
                    double fitT0_0 = rawhit->getT0(0);
                    double fitT0_1 = rawhit->getT0(1);
                    int fitN = rawhit->getFitN();
                    double doublepulseT0 = 0.0;
                    int inTimeFit = -1;

                    int cluster_count = hit->getRawHits()->GetEntries();

                    if(fitN == 1)
                    {
                        double fitAmp = rawhit->getAmp(0);
                        trkHistos_->Fill1DHisto(charge_str+"_trk_singlepulse_t0_h", fitT0_0, weight);
                        trkHistos_->Fill2DHisto(charge_str+"_trk_singlepulse_t0_v_layer_hh", layer, fitT0_0, weight);
                        trkHistos_->Fill2DHisto(charge_str+"_trk_singlepulse_amplitude_v_t0_hh",fitT0_0,fitAmp, weight);

                        trkHistos_->Fill2DHisto(charge_str+"_trk_"+std::to_string(cluster_count)+"_hitcluster_rawhit_amplitude_v_t0_hh",fitT0_0,fitAmp, weight);
                    }

                    else if(fitN == 2){
                        if(std::abs(0-fitT0_0) < std::abs(0-fitT0_1)){
                            inTimeFit = 0;
                            doublepulseT0 = fitT0_0;
                        }
                        else{
                            inTimeFit = 1;
                            doublepulseT0 = fitT0_1;
                        }
                        double fitAmp = rawhit->getAmp(inTimeFit);

                        trkHistos_->Fill1DHisto(charge_str+"_trk_doublepulse_t0_h", doublepulseT0, weight);
                        trkHistos_->Fill2DHisto(charge_str+"_trk_doublepulse_t0_v_layer_hh", layer, doublepulseT0, weight);
                        trkHistos_->Fill2DHisto(charge_str+"_trk_doublepulse_amplitude_v_t0_hh", doublepulseT0, fitAmp, weight);
                        trkHistos_->Fill2DHisto(charge_str+"_trk_"+std::to_string(cluster_count)+"_hitcluster_rawhit_amplitude_v_t0_hh",fitT0_0,fitAmp, weight);
                    }
                }
            }

        }
        
        trkHistos_->Fill1DHistograms(track);
        trkHistos_->Fill2DTrack(track);

        /*
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
        */
        
        n_sel_tracks++;
    }//Loop on tracks

    //trkHistos_->Fill1DHisto("n_tracks_h",n_sel_tracks);

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
