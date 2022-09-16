/**
 *@file VertexAnaProcessor.cxx
 *@brief Main vertex analysis processor
 *@author PF, SLAC
 */

#include "VertexAnaProcessor.h"
#include <iostream>
#include <map>
#include <math.h> 
#include "TRandom.h"

VertexAnaProcessor::VertexAnaProcessor(const std::string& name, Process& process) : Processor(name,process) {

}

//TODO Check this destructor

VertexAnaProcessor::~VertexAnaProcessor(){}

void VertexAnaProcessor::configure(const ParameterSet& parameters) {
    std::cout << "Configuring VertexAnaProcessor" <<std::endl;
    try
    {
        debug_   = parameters.getInteger("debug",debug_);
        anaName_ = parameters.getString("anaName",anaName_);
        tsColl_  = parameters.getString("tsColl",tsColl_);
        vtxColl_ = parameters.getString("vtxColl",vtxColl_);
        trkColl_ = parameters.getString("trkColl",trkColl_);
        hitColl_ = parameters.getString("hitColl",hitColl_);
        ecalColl_ = parameters.getString("ecalColl",ecalColl_);
        mcColl_  = parameters.getString("mcColl",mcColl_);
        isRadPDG_ = parameters.getInteger("isRadPDG",isRadPDG_);
        smearMCTrackTime_ = parameters.getInteger("smearMCTrackTime",smearMCTrackTime_);
        std::cout << "smearMCTrackTime: " << smearMCTrackTime_ << std::endl;

        selectionCfg_   = parameters.getString("vtxSelectionjson",selectionCfg_);
        histoCfg_ = parameters.getString("histoCfg",histoCfg_);
        mcHistoCfg_ = parameters.getString("mcHistoCfg",mcHistoCfg_);
        timeOffset_ = parameters.getDouble("CalTimeOffset",timeOffset_);
        beamE_  = parameters.getDouble("beamE",beamE_);
        isData_  = parameters.getInteger("isData",isData_);
        analysis_        = parameters.getString("analysis");

        //region definitions
        regionSelections_ = parameters.getVString("regionDefinitions",regionSelections_);


    }
    catch (std::runtime_error& error)
    {
        std::cout<<error.what()<<std::endl;
    }
}

void VertexAnaProcessor::initialize(TTree* tree) {
    tree_ = tree;
    _ah =  std::make_shared<AnaHelpers>();

    //vtxSelector  = std::make_shared<BaseSelector>(anaName_+"_"+"vtxSelection",selectionCfg_);
    //vtxSelector->setDebug(debug_);
    //vtxSelector->LoadSelection();

    //_vtx_histos = std::make_shared<TrackHistos>(anaName_+"_"+"vtxSelection");
    //_vtx_histos->loadHistoConfig(histoCfg_);
    //_vtx_histos->DefineHistos();

    if(!isData_){
        _mc_vtx_histos = std::make_shared<MCAnaHistos>(anaName_+"_mc_"+"vtxSelection");
        _mc_vtx_histos->loadHistoConfig(mcHistoCfg_);
        _mc_vtx_histos->DefineHistos();
        _mc_vtx_histos->Define2DHistos();
    }


    //    histos = new MCAnaHistos(anaName_);
    //histos->loadHistoConfig(histCfgFilename_)
    //histos->DefineHistos();
    //histos->Define2DHistos();


    //For each region initialize plots

    for (unsigned int i_reg = 0; i_reg < regionSelections_.size(); i_reg++) {
        std::string regname = AnaHelpers::getFileName(regionSelections_[i_reg],false);
        std::cout<<"Setting up region:: " << regname <<std::endl;
        _reg_vtx_selectors[regname] = std::make_shared<BaseSelector>(anaName_+"_"+regname, regionSelections_[i_reg]);
        _reg_vtx_selectors[regname]->setDebug(debug_);
        _reg_vtx_selectors[regname]->LoadSelection();

        _reg_vtx_histos[regname] = std::make_shared<TrackHistos>(anaName_+"_"+regname);
        _reg_vtx_histos[regname]->loadHistoConfig(histoCfg_);
        _reg_vtx_histos[regname]->DefineHistos();


        _reg_mc_vtx_histos[regname] = std::make_shared<MCAnaHistos>(anaName_+"_mc_"+regname);
        _reg_mc_vtx_histos[regname]->loadHistoConfig(mcHistoCfg_);
        _reg_mc_vtx_histos[regname]->DefineHistos();



        _reg_tuples[regname] = std::make_shared<FlatTupleMaker>(anaName_+"_"+regname+"_tree");
        _reg_tuples[regname]->addVariable("unc_vtx_mass");
        _reg_tuples[regname]->addVariable("unc_vtx_z");
        if(!isData_)
        {
            _reg_tuples[regname]->addVariable("true_vtx_z");
            _reg_tuples[regname]->addVariable("true_vtx_mass");
        }

        _regions.push_back(regname);
    }

    // Get list of branches in tree to help protect accessing them
    int nBr = tree_->GetListOfBranches()->GetEntries();
    if (debug_) std::cout << "Tree has " << nBr << " branches" << std::endl;
    for(int iBr = 0; iBr < nBr; iBr++)
    {
        TBranch *br = dynamic_cast<TBranch*>(tree_->GetListOfBranches()->At(iBr));
        brMap_.insert(std::map<const char *, int, char_cmp>::value_type(br->GetName(), 1));
        if (debug_) std::cout << br->GetName() << ": " << brMap_[br->GetName()] << std::endl;
    }

    //init Reading Tree
    tree_->SetBranchAddress("EventHeader", &evth_ , &bevth_);
    if (brMap_.find(tsColl_.c_str()) != brMap_.end()) tree_->SetBranchAddress(tsColl_.c_str(), &ts_ , &bts_);
    tree_->SetBranchAddress(vtxColl_.c_str(), &vtxs_ , &bvtxs_);
    tree_->SetBranchAddress(hitColl_.c_str(), &hits_   , &bhits_);
    tree_->SetBranchAddress(ecalColl_.c_str(), &ecal_  , &becal_);
    if(!isData_ && !mcColl_.empty()) tree_->SetBranchAddress(mcColl_.c_str() , &mcParts_, &bmcParts_);
    //If track collection name is empty take the tracks from the particles. TODO:: change this
    if (!trkColl_.empty())
        tree_->SetBranchAddress(trkColl_.c_str(),&trks_, &btrks_);
}

bool VertexAnaProcessor::process(IEvent* ievent) {
    if(debug_) {
        std:: cout << "----------------- Event " << evth_->getEventNumber() << " -----------------" << std::endl;
    }
    HpsEvent* hps_evt = (HpsEvent*) ievent;
    double weight = 1.;


    //Get "true" mass
    double apMass = -0.9;
    double apZ = -0.9;

    //Plot info about which trigger bits are present in the event
    if (ts_ != nullptr)
    {
        //_vtx_histos->Fill2DHisto("trig_count_hh", 
        //        ((int)ts_->prescaled.Single_3_Top)+((int)ts_->prescaled.Single_3_Bot),
        //        ((int)ts_->prescaled.Single_2_Top)+((int)ts_->prescaled.Single_2_Bot));
    }
    int NposTrks = 0;
    int NeleTrks = 0;
    for (int iT = 0; iT < trks_->size(); iT++)
    {
        if (trks_->at(iT)->getCharge() > 0) NposTrks++;
        else NeleTrks++;
    }
    //_vtx_histos->Fill2DHisto("n_tracks_hh", NeleTrks, NposTrks); 
    //_vtx_histos->Fill1DHisto("n_vtx_h", vtxs_->size()); 

    if (mcParts_) {
        for(int i = 0; i < mcParts_->size(); i++)
        {
            if(mcParts_->at(i)->getPDG() == 622)
            {
                apMass = mcParts_->at(i)->getMass();
                apZ = mcParts_->at(i)->getVertexPosition().at(2);
            }
        }

        if (!isData_) _mc_vtx_histos->FillMCParticles(mcParts_, analysis_);
    }
    //Store processed number of events
    std::vector<Vertex*> selected_vtxs;
    bool passVtxPresel = false;

    // Fill some diagnostic histos
    /*
    for ( int i_ecal = 0; i_ecal < ecal_->size(); i_ecal++ ) {

        if (vtxs_->size() == 0){
            _vtx_histos->Fill1DHisto("EecalClus_noVtxs_h",ecal_->at(i_ecal)->getEnergy());
        } else {
            _vtx_histos->Fill1DHisto("EecalClus_isVtxs_h",ecal_->at(i_ecal)->getEnergy());
        }
    }


    if (vtxs_->size() == 0){
        _vtx_histos->Fill1DHisto("n_ecalClus_noVtxs_h",ecal_->size());
        _vtx_histos->Fill1DHisto("n_tracks_noVtxs_h",trks_->size());
        for (int i_trk = 0; i_trk < trks_->size(); i_trk++ ){
            _vtx_histos->Fill1DHisto("Ptracks_noVtxs_h",trks_->at(i_trk)->getP());
        }

    } else {
        _vtx_histos->Fill1DHisto("n_ecalClus_isVtxs_h",ecal_->size());
        _vtx_histos->Fill1DHisto("n_tracks_isVtxs_h",trks_->size());
        for (int i_trk = 0; i_trk < trks_->size(); i_trk++ ){
            _vtx_histos->Fill1DHisto("Ptracks_isVtxs_h",trks_->at(i_trk)->getP());
        }
    }

    if(debug_){
        std::cout<<"Number of vertices found in event: "<< vtxs_->size()<<std::endl;
    }

    _vtx_histos->Fill1DHisto("n_vertices_h",selected_vtxs.size());
    if (trks_)
        _vtx_histos->Fill1DHisto("n_tracks_h",trks_->size()); 
        */

    //Make Plots for each region: loop on each region. Check if the region has the cut and apply it
    //TODO Clean this up => Cuts should be implemented in each region?
    //TODO Bring the preselection out of this stupid loop

    //TODO add yields. => Quite terrible way to loop.
    for (auto region : _regions ) {

        int nGoodVtx = 0;
        Vertex* goodVtx = nullptr;

        float truePsum = -1;
        float trueEsum = -1;

        for ( int i_vtx = 0; i_vtx <  vtxs_->size(); i_vtx++ ) {

            Vertex* vtx = vtxs_->at(i_vtx);

            //No cuts.
            _reg_vtx_selectors[region]->getCutFlowHisto()->Fill(0.,weight);

            Particle* ele = nullptr;
            Particle* pos = nullptr;

            _ah->GetParticlesFromVtx(vtx,ele,pos);

            CalCluster eleClus = ele->getCluster();
            CalCluster posClus = pos->getCluster();

            double ele_E = ele->getEnergy();
            double pos_E = pos->getEnergy();

            //Compute analysis variables here.

            Track ele_trk = ele->getTrack();
            Track pos_trk = pos->getTrack();

            //Get the shared info - TODO change and improve

            Track* ele_trk_gbl = nullptr;
            Track* pos_trk_gbl = nullptr;

            if (!trkColl_.empty()) {
                bool foundTracks = _ah->MatchToGBLTracks(ele_trk.getID(),pos_trk.getID(),
                        ele_trk_gbl, pos_trk_gbl, *trks_);

                if (!foundTracks) {
                    if (debug_)
                        std::cout<<"VertexAnaProcessor::ERROR couldn't find ele/pos in the "<<trkColl_ <<"collection"<<std::endl;
                    continue;
                }
            }
            else {

                ele_trk_gbl = (Track*) ele_trk.Clone();
                pos_trk_gbl = (Track*) pos_trk.Clone();
            }

            //smear MC Track time
            if(!isData_ && smearMCTrackTime_){
                smearMCTrackTime(ele_trk_gbl);
                smearMCTrackTime(pos_trk_gbl);
            }

            //Add the momenta to the tracks
            //ele_trk_gbl->setMomentum(ele->getMomentum()[0],ele->getMomentum()[1],ele->getMomentum()[2]);
            //pos_trk_gbl->setMomentum(pos->getMomentum()[0],pos->getMomentum()[1],pos->getMomentum()[2]);
            TVector3 recEleP(ele->getMomentum()[0],ele->getMomentum()[1],ele->getMomentum()[2]);
            TLorentzVector p_ele;
            p_ele.SetPxPyPzE(ele_trk_gbl->getMomentum()[0],ele_trk_gbl->getMomentum()[1],ele_trk_gbl->getMomentum()[2], ele_E);
            TLorentzVector p_pos;
            p_pos.SetPxPyPzE(pos_trk_gbl->getMomentum()[0],pos_trk_gbl->getMomentum()[1],pos_trk_gbl->getMomentum()[2], pos_E);

            //Defining these here so they are in scope elsewhere
            TVector3 trueEleP;
            TVector3 truePosP;

            if (debug_) {
                std::cout<<"Check on ele_Track"<<std::endl;
                std::cout<<"Number of hits:"<<ele_trk_gbl->getTrackerHitCount()<<std::endl;
            }

            bool foundL1ele = false;
            bool foundL2ele = false;
            _ah->InnermostLayerCheck(ele_trk_gbl, foundL1ele, foundL2ele);


            if (debug_) {
                std::cout<<"Check on pos_Track"<<std::endl;
                std::cout<<"Number of hits:"<<ele_trk_gbl->getTrackerHitCount()<<std::endl;
            }
            bool foundL1pos = false;
            bool foundL2pos = false;

            _ah->InnermostLayerCheck(pos_trk_gbl, foundL1pos, foundL2pos);

            if (debug_) {
                std::cout<<"Check on pos_Track"<<std::endl;
                std::cout<<"Innermost:"<<foundL1pos<<" Second Innermost:"<<foundL2pos<<std::endl;
            }

            double corr_eleClusterTime = ele->getCluster().getTime() - timeOffset_;
            double corr_posClusterTime = pos->getCluster().getTime() - timeOffset_;

            double botClusTime = 0.0;
            if(ele->getCluster().getPosition().at(1) < 0.0) botClusTime = ele->getCluster().getTime();
            else botClusTime = pos->getCluster().getTime();

            //Momentum
            TVector3 ele_mom;
            ele_mom.SetX(ele_trk_gbl->getMomentum()[0]);
            ele_mom.SetY(ele_trk_gbl->getMomentum()[1]);
            ele_mom.SetZ(ele_trk_gbl->getMomentum()[2]);


            TVector3 pos_mom;
            pos_mom.SetX(pos_trk_gbl->getMomentum()[0]);
            pos_mom.SetY(pos_trk_gbl->getMomentum()[1]);
            pos_mom.SetZ(pos_trk_gbl->getMomentum()[2]);
            
            //Pair 1 cut
            if (isData_) {
                if (!_reg_vtx_selectors[region]->passCutEq("Pair1_eq",(int)evth_->isPair1Trigger(),weight))
                    break;
            }

            if (!_reg_vtx_selectors[region]->passCutLt("eleTrkTime_lt",fabs(ele_trk_gbl->getTrackTime()),weight))
                continue;

            if (!_reg_vtx_selectors[region]->passCutLt("posTrkTime_lt",fabs(pos_trk_gbl->getTrackTime()),weight))
                continue;

            //Less than 4 shared hits for ele/pos track
            if (!_reg_vtx_selectors[region]->passCutLt("eleNshared_lt",ele_trk_gbl->getNShared(),weight)) {
                continue;
            }

            if (!_reg_vtx_selectors[region]->passCutLt("posNshared_lt",pos_trk_gbl->getNShared(),weight)) {
                continue;
            }

            if (!_reg_vtx_selectors[region]->passCutLt("eleTrkCluMatch_lt",ele->getGoodnessOfPID(),weight))
                continue;

            //Pos Track-cluster match
            if (!_reg_vtx_selectors[region]->passCutLt("posTrkCluMatch_lt",pos->getGoodnessOfPID(),weight))
                continue;

            //Ele Pos Cluster Time Difference
            if (!_reg_vtx_selectors[region]->passCutLt("eleposCluTimeDiff_lt",fabs(corr_eleClusterTime - corr_posClusterTime),weight))
                continue;

            //Ele Track-Cluster Time Difference
            if (!_reg_vtx_selectors[region]->passCutLt("eleTrkCluTimeDiff_lt",fabs(ele_trk_gbl->getTrackTime() - corr_eleClusterTime),weight))
                continue;

            //Pos Track-Cluster Time Difference
            if (!_reg_vtx_selectors[region]->passCutLt("posTrkCluTimeDiff_lt",fabs(pos_trk_gbl->getTrackTime() - corr_posClusterTime),weight))
                continue;

            //Beam Electron cut
            if (!_reg_vtx_selectors[region]->passCutLt("eleMom_lt",ele_mom.Mag(),weight))
                continue;

            //PSum low cut
            if (!_reg_vtx_selectors[region]->passCutLt("pSum_lt",(p_ele.P()+p_pos.P()),weight))
                continue;

            //Ele Track Quality - Chi2
            if (!_reg_vtx_selectors[region]->passCutLt("eleTrkChi2_lt",ele_trk_gbl->getChi2(),weight))
                continue;

            //Pos Track Quality - Chi2
            if (!_reg_vtx_selectors[region]->passCutLt("posTrkChi2_lt",pos_trk_gbl->getChi2(),weight))
                continue;

            //Ele Track Quality - Chi2Ndf
            if (!_reg_vtx_selectors[region]->passCutLt("eleTrkChi2Ndf_lt",ele_trk_gbl->getChi2Ndf(),weight))
                continue;

            //Pos Track Quality - Chi2Ndf
            if (!_reg_vtx_selectors[region]->passCutLt("posTrkChi2Ndf_lt",pos_trk_gbl->getChi2Ndf(),weight))
                continue;


            //Ele min momentum cut
            if (!_reg_vtx_selectors[region]->passCutGt("eleMom_gt",ele_mom.Mag(),weight))
                continue;

            //Pos min momentum cut
            if (!_reg_vtx_selectors[region]->passCutGt("posMom_gt",pos_mom.Mag(),weight))
                continue;

            //Ele nHits
            int ele2dHits = ele_trk_gbl->getTrackerHitCount();
            if (!ele_trk_gbl->isKalmanTrack())
                ele2dHits*=2;

            if (!_reg_vtx_selectors[region]->passCutGt("eleN2Dhits_gt",ele2dHits,weight))  {
                continue;
            }

            //Pos nHits
             int pos2dHits = pos_trk_gbl->getTrackerHitCount();
             if (!pos_trk_gbl->isKalmanTrack())
                pos2dHits*=2;

             if (!_reg_vtx_selectors[region]->passCutGt("posN2Dhits_gt",pos2dHits,weight))  {
                continue;
             }

            //Min vtx momentum
            if (!_reg_vtx_selectors[region]->passCutGt("minVtxMom_gt",(ele_mom+pos_mom).Mag(),weight))
                continue;

            //vtx Z position
            if (!_reg_vtx_selectors[region]->passCutGt("uncVtxZ_gt",vtx->getZ(),weight))
                continue;

            //Chi2
            if (!_reg_vtx_selectors[region]->passCutLt("chi2unc_lt",vtx->getChi2(),weight))
                continue;

            //PSum high cut
            if (!_reg_vtx_selectors[region]->passCutGt("pSum_gt",(p_ele.P()+p_pos.P()),weight))
                continue;

            //L1 requirement
            if (!_reg_vtx_selectors[region]->passCutEq("L1Requirement_eq",(int)(foundL1ele&&foundL1pos),weight))
                continue;

            //L2 requirement
            if (!_reg_vtx_selectors[region]->passCutEq("L2Requirement_eq",(int)(foundL2ele&&foundL2pos),weight))
                continue;

            //L1 requirement for positron
            if (!_reg_vtx_selectors[region]->passCutEq("L1PosReq_eq",(int)(foundL1pos),weight))
                continue;

            //Require Electron Cluster exists
            if (!_reg_vtx_selectors[region]->passCutGt("eleClusE_gt",eleClus.getEnergy(),weight))
                continue;

            //Require Electron Cluster does NOT exists
            if (!_reg_vtx_selectors[region]->passCutLt("eleClusE_lt",eleClus.getEnergy(),weight))
                continue;

            //Require Positron Cluster exists
            if (!_reg_vtx_selectors[region]->passCutGt("posClusE_gt",posClus.getEnergy(),weight))
                continue;

            //Require Positron Cluster does NOT exists
            if (!_reg_vtx_selectors[region]->passCutLt("posClusE_lt",posClus.getEnergy(),weight))
                continue;

            //Tracking Volume for positron
            if (!_reg_vtx_selectors[region]->passCutGt("volPos_top", p_pos.Py(), weight))
                continue;

            if (!_reg_vtx_selectors[region]->passCutLt("volPos_bot", p_pos.Py(), weight))
                continue;


            //Bottom Cluster Time
            if (!_reg_vtx_selectors[region]->passCutLt("botCluTime_lt", botClusTime, weight))
                continue;

            if (!_reg_vtx_selectors[region]->passCutGt("botCluTime_gt", botClusTime, weight))
                continue;


            //ESum low cut
            if (!_reg_vtx_selectors[region]->passCutLt("eSum_lt",(ele_E+pos_E),weight))
                continue;

            //ESum high cut
            if (!_reg_vtx_selectors[region]->passCutGt("eSum_gt",(ele_E+pos_E),weight))
                continue;

            //No shared hits requirement
            if (!_reg_vtx_selectors[region]->passCutEq("ele_sharedL0_eq",(int)ele_trk_gbl->getSharedLy0(),weight))
                continue;
            if (!_reg_vtx_selectors[region]->passCutEq("pos_sharedL0_eq",(int)pos_trk_gbl->getSharedLy0(),weight))
                continue;
            if (!_reg_vtx_selectors[region]->passCutEq("ele_sharedL1_eq",(int)ele_trk_gbl->getSharedLy1(),weight))
                continue;
            if (!_reg_vtx_selectors[region]->passCutEq("pos_sharedL1_eq",(int)pos_trk_gbl->getSharedLy1(),weight))
                continue;

            //Min vtx Y pos
            if (!_reg_vtx_selectors[region]->passCutGt("VtxYPos_gt", vtx->getY(), weight))
                continue;

            //Max vtx Y pos
            if (!_reg_vtx_selectors[region]->passCutLt("VtxYPos_lt", vtx->getY(), weight))
                continue;

            //If this is MC check if MCParticle matched to the electron track is from rad or recoil
            if(!isData_)
            {

                //Fill MC plots after all selections
                if (!isData_) _reg_mc_vtx_histos[region]->FillMCParticles(mcParts_, analysis_);

                //Build map of hits and the associated MC part ids for later
                TRefArray* ele_trk_hits = ele_trk_gbl->getSvtHits();
                std::map<int, std::vector<int> > trueHitIDs;
                for(int i = 0; i < hits_->size(); i++)
                {
                    TrackerHit* hit = hits_->at(i);
                    trueHitIDs[hit->getID()] = hit->getMCPartIDs();
                }
                //std::cout << "There are " << ele_trk_hits->GetEntries() << " hits on this track" << std::endl;
                //Count the number of hits per part on the track
                std::map<int, int> nHits4part;
                for(int i = 0; i < ele_trk_hits->GetEntries(); i++)
                {
                    TrackerHit* eleHit = (TrackerHit*)ele_trk_hits->At(i);
                    for(int idI = 0; idI < trueHitIDs[eleHit->getID()].size(); idI++ )
                    {
                        int partID = trueHitIDs[eleHit->getID()].at(idI);
                        if ( nHits4part.find(partID) == nHits4part.end() )
                        {
                            // not found
                            nHits4part[partID] = 1;
                        }
                        else
                        {
                            // found
                            nHits4part[partID]++;
                        }
                    }
                }

                //Determine the MC part with the most hits on the track
                int maxNHits = 0;
                int maxID = 0;
                for (std::map<int,int>::iterator it=nHits4part.begin(); it!=nHits4part.end(); ++it)
                {
                    if(it->second > maxNHits)
                    {
                        maxNHits = it->second;
                        maxID = it->first;
                    }
                }

                //Find the correct mc part and grab mother id
                int isRadEle = -999;
                int isRecEle = -999;


                trueEleP.SetXYZ(-999,-999,-999);
                truePosP.SetXYZ(-999,-999,-999);
                if (mcParts_) {
                    float trueEleE = -1;
                    float truePosE = -1;
                    for(int i = 0; i < mcParts_->size(); i++)
                    {
                        int momPDG = mcParts_->at(i)->getMomPDG();
                        if(mcParts_->at(i)->getPDG() == 11 && momPDG == 622)
                        {
                            std::vector<double> lP = mcParts_->at(i)->getMomentum();
                            trueEleP.SetXYZ(lP[0],lP[1],lP[2]);
                            trueEleE = mcParts_->at(i)->getEnergy();

                        }
                        if(mcParts_->at(i)->getPDG() == -11 && momPDG == 622)
                        {
                            std::vector<double> lP = mcParts_->at(i)->getMomentum();
                            truePosP.SetXYZ(lP[0],lP[1],lP[2]);
                            truePosE = mcParts_->at(i)->getEnergy();

                        }
                        if(trueEleP.X() != -999 && truePosP.X() != -999){
                            truePsum =  trueEleP.Mag() + trueEleP.Mag();
                            trueEsum = trueEleE + truePosE;
                        }

                        if(mcParts_->at(i)->getID() != maxID) continue;
                        //Default isRadPDG = 622
                        if(momPDG == isRadPDG_) isRadEle = 1;
                        if(momPDG == 623) isRecEle = 1;
                    }
                }
                double momRatio = recEleP.Mag() / trueEleP.Mag();
                double momAngle = trueEleP.Angle(recEleP) * TMath::RadToDeg();
                if (!_reg_vtx_selectors[region]->passCutLt("momRatio_lt", momRatio, weight)) continue;
                if (!_reg_vtx_selectors[region]->passCutGt("momRatio_gt", momRatio, weight)) continue;
                if (!_reg_vtx_selectors[region]->passCutLt("momAngle_lt", momAngle, weight)) continue;

                if (!_reg_vtx_selectors[region]->passCutEq("isRadEle_eq", isRadEle, weight)) continue;
                if (!_reg_vtx_selectors[region]->passCutEq("isRecEle_eq", isRecEle, weight)) continue;
            }

            goodVtx = vtx;
            nGoodVtx++;
        }

        //N selected vertices - this is quite a silly cut to make at the end. But okay. that's how we decided atm.
        if (!_reg_vtx_selectors[region]->passCutEq("nVtxs_eq", nGoodVtx, weight))
            continue;
        //Move to after N vertices cut (was filled before)
        _reg_vtx_histos[region]->Fill1DHisto("n_vertices_h", nGoodVtx, weight);

        Vertex* vtx = goodVtx;

        Particle* ele = nullptr;
        Particle* pos = nullptr;

        if (!vtx || !_ah->GetParticlesFromVtx(vtx,ele,pos))
            continue;

        CalCluster eleClus = ele->getCluster();
        CalCluster posClus = pos->getCluster();

        double ele_E = ele->getEnergy();
        double pos_E = pos->getEnergy();

        //Compute analysis variables here.
        Track ele_trk = ele->getTrack();
        Track pos_trk = pos->getTrack();
        //Get the shared info - TODO change and improve

        Track* ele_trk_gbl = nullptr;
        Track* pos_trk_gbl = nullptr;

        if (!trkColl_.empty()) {
            bool foundTracks = _ah->MatchToGBLTracks(ele_trk.getID(),pos_trk.getID(),
                    ele_trk_gbl, pos_trk_gbl, *trks_);

            if (!foundTracks) {
                if (debug_)
                    std::cout<<"VertexAnaProcessor::ERROR couldn't find ele/pos in the "<<trkColl_ <<"collection"<<std::endl;
                continue;
            }
        }
        else {

            ele_trk_gbl = (Track*) ele_trk.Clone();
            pos_trk_gbl = (Track*) pos_trk.Clone();
        }
        if(ts_ != nullptr)
        {
            _reg_vtx_histos[region]->Fill2DHisto("trig_count_hh", 
                    ((int)ts_->prescaled.Single_3_Top)+((int)ts_->prescaled.Single_3_Bot),
                    ((int)ts_->prescaled.Single_2_Top)+((int)ts_->prescaled.Single_2_Bot));
        }
        _reg_vtx_histos[region]->Fill1DHisto("n_vtx_h", vtxs_->size()); 
        int NposTrks = 0;
        int NeleTrks = 0;
        for (int iT = 0; iT < trks_->size(); iT++)
        {
            if (trks_->at(iT)->getCharge() > 0) NposTrks++;
            else NeleTrks++;
        }

        TVector3 recEleP(ele->getMomentum()[0],ele->getMomentum()[1],ele->getMomentum()[2]);
        TLorentzVector p_ele;
        p_ele.SetPxPyPzE(ele_trk_gbl->getMomentum()[0],ele_trk_gbl->getMomentum()[1],ele_trk_gbl->getMomentum()[2], ele_E);
        TLorentzVector p_pos;
        p_pos.SetPxPyPzE(pos_trk_gbl->getMomentum()[0],pos_trk_gbl->getMomentum()[1],pos_trk_gbl->getMomentum()[2], pos_E);

        //Trackhits on Layers
        std::vector<int> hit_layers;
        int hitCode = 0;
        for (int ihit = 0; ihit<ele_trk_gbl->getSvtHits()->GetEntries(); ++ihit) {
            TrackerHit* hit = (TrackerHit*) ele_trk_gbl->getSvtHits()->At(ihit);
            int layer = hit->getLayer();
            _reg_vtx_histos[region]->Fill1DHisto("ele_hitlayers_h", layer);
            _reg_vtx_histos[region]->Fill2DHisto("ele_tanlambda_v_hitlayers_hh",layer,ele_trk_gbl->getTanLambda(), weight);
            _reg_vtx_histos[region]->Fill2DHisto("ele_phi_v_hitlayers_hh",layer,ele_trk_gbl->getPhi(), weight);
            _reg_vtx_histos[region]->Fill2DHisto("ele_p_v_hitlayers_hh",layer,p_ele.P(), weight);
        }
        for (int ihit = 0; ihit<pos_trk_gbl->getSvtHits()->GetEntries(); ++ihit) {
            TrackerHit* hit = (TrackerHit*) pos_trk_gbl->getSvtHits()->At(ihit);
            int layer = hit->getLayer();
            _reg_vtx_histos[region]->Fill1DHisto("pos_hitlayers_h", layer);
            _reg_vtx_histos[region]->Fill2DHisto("pos_tanlambda_v_hitlayers_hh",layer,pos_trk_gbl->getTanLambda(), weight);
            _reg_vtx_histos[region]->Fill2DHisto("pos_phi_v_hitlayers_hh",layer,pos_trk_gbl->getPhi(), weight);
            _reg_vtx_histos[region]->Fill2DHisto("pos_p_v_hitlayers_hh",layer,p_pos.P(), weight);
        }

        _reg_vtx_histos[region]->Fill2DHisto("n_tracks_hh", NeleTrks, NposTrks); 

        //Add the momenta to the tracks
        //ele_trk_gbl->setMomentum(ele->getMomentum()[0],ele->getMomentum()[1],ele->getMomentum()[2]);
        //pos_trk_gbl->setMomentum(pos->getMomentum()[0],pos->getMomentum()[1],pos->getMomentum()[2]);
        //

        double corr_eleClusterTime = ele->getCluster().getTime() - timeOffset_;
        double corr_posClusterTime = pos->getCluster().getTime() - timeOffset_;

        //Post selection plots
        
        //Track Cluster matcher plots
        _reg_vtx_histos[region]->Fill1DHisto("ele_track_clus_dt_h",ele_trk_gbl->getTrackTime() - corr_eleClusterTime, weight);
        _reg_vtx_histos[region]->Fill1DHisto("pos_track_clus_dt_h",pos_trk_gbl->getTrackTime() - corr_posClusterTime, weight);
        _reg_vtx_histos[region]->Fill1DHisto("ele_pos_clusTimeDiff_h",corr_eleClusterTime - corr_posClusterTime),weight;
         //2d histos
        _reg_vtx_histos[region]->Fill2DHisto("ele_track_clus_dt_v_p_hh",p_ele.P(),ele_trk_gbl->getTrackTime() - corr_eleClusterTime, weight);
        _reg_vtx_histos[region]->Fill2DHisto("pos_track_clus_dt_v_p_hh",p_pos.P(),pos_trk_gbl->getTrackTime() - corr_posClusterTime, weight);
        _reg_vtx_histos[region]->Fill2DHisto("ele_pos_clusTimeDiff_v_pSum_hh",p_ele.P()+p_pos.P(),corr_eleClusterTime - corr_posClusterTime,weight);
        _reg_vtx_histos[region]->Fill2DHisto("ele_clusT_v_ele_trackT_hh",ele_trk_gbl->getTrackTime(),corr_eleClusterTime, weight);
        _reg_vtx_histos[region]->Fill2DHisto("pos_clusT_v_pos_trackT_hh",pos_trk_gbl->getTrackTime(),corr_posClusterTime, weight);
        _reg_vtx_histos[region]->Fill2DHisto("pos_v_ele_track_p_hh",p_ele.P(),p_pos.P(), weight);

        _reg_vtx_histos[region]->Fill2DHistograms(vtx,weight);
        _reg_vtx_histos[region]->Fill1DVertex(vtx,
                ele,
                pos,
                ele_trk_gbl,
                pos_trk_gbl,
                weight);

        _reg_vtx_histos[region]->Fill1DHisto("vtx_Psum_h", p_ele.P()+p_pos.P(), weight);
        _reg_vtx_histos[region]->Fill1DHisto("vtx_Esum_h", eleClus.getEnergy()+posClus.getEnergy(), weight);
        _reg_vtx_histos[region]->Fill2DHisto("ele_vtxZ_iso_hh", TMath::Min(ele_trk_gbl->getIsolation(0), ele_trk_gbl->getIsolation(1)), vtx->getZ(), weight);
        _reg_vtx_histos[region]->Fill2DHisto("pos_vtxZ_iso_hh", TMath::Min(pos_trk_gbl->getIsolation(0), pos_trk_gbl->getIsolation(1)), vtx->getZ(), weight);
        _reg_vtx_histos[region]->Fill2DTrack(ele_trk_gbl,weight,"ele_");
        _reg_vtx_histos[region]->Fill2DTrack(pos_trk_gbl,weight,"pos_");
        _reg_vtx_histos[region]->Fill1DHisto("mcMass622_h",apMass);
        _reg_vtx_histos[region]->Fill1DHisto("mcZ622_h",apZ);

        if (trks_) _reg_vtx_histos[region]->Fill1DHisto("n_tracks_h",trks_->size(),weight);

        //Just for the selected vertex
        _reg_tuples[region]->setVariableValue("unc_vtx_mass", vtx->getInvMass());
        if(!isData_)
        {
            _reg_vtx_histos[region]->Fill2DHisto("vtx_Esum_vs_true_Esum_hh",eleClus.getEnergy()+posClus.getEnergy(), trueEsum, weight);
            _reg_vtx_histos[region]->Fill2DHisto("vtx_Psum_vs_true_Psum_hh",p_ele.P()+p_pos.P(), truePsum, weight);
            _reg_tuples[region]->setVariableValue("true_vtx_z", apZ);
            _reg_tuples[region]->setVariableValue("true_vtx_mass", apMass);
            _reg_vtx_histos[region]->Fill1DHisto("true_vtx_psum_h",truePsum,weight);
        }

        //TODO put this in the Vertex!
        TVector3 vtxPosSvt;
        vtxPosSvt.SetX(vtx->getX());
        vtxPosSvt.SetY(vtx->getY());
        vtxPosSvt.SetZ(vtx->getZ());
        vtxPosSvt.RotateY(-0.0305);

        _reg_tuples[region]->setVariableValue("unc_vtx_z"   , vtxPosSvt.Z());
        _reg_tuples[region]->fill();
    }// regions



    return true;
}

void VertexAnaProcessor::smearMCTrackTime(Track* track){

    //smear track time in MC
    double mc_tres = 0.62;
    double data_tres = 2.1;
    double smear_tres = sqrt( (data_tres*data_tres) - (mc_tres*mc_tres) );
    TRandom* rand = new TRandom();
    rand->SetSeed();
    double rand_smear = rand->Gaus(0,smear_tres);
    double track_time = track->getTrackTime();
    track->setTrackTime(track_time + rand_smear);
}

void VertexAnaProcessor::finalize() {

    //TODO clean this up a little.
    outF_->cd();
    //_vtx_histos->saveHistos(outF_,_vtx_histos->getName());
    //outF_->cd(_vtx_histos->getName().c_str());
    //vtxSelector->getCutFlowHisto()->Write();

    outF_->cd();
    if(!isData_)
        _mc_vtx_histos->saveHistos(outF_, _mc_vtx_histos->getName());
    //delete histos;
    //histos = nullptr;


    for (reg_it it = _reg_vtx_histos.begin(); it!=_reg_vtx_histos.end(); ++it) {
        std::string dirName = anaName_+"_"+it->first;
        (it->second)->saveHistos(outF_,dirName);
        outF_->cd(dirName.c_str());
        _reg_vtx_selectors[it->first]->getCutFlowHisto()->Write();
        //Save tuples
        _reg_tuples[it->first]->writeTree();
    }

    for (reg_mc_it it = _reg_mc_vtx_histos.begin(); it!=_reg_mc_vtx_histos.end(); ++it) {
        std::string dirName = anaName_+"_mc_"+it->first;
        (it->second)->saveHistos(outF_,dirName);
        outF_->cd(dirName.c_str());
    }

    outF_->Close();

}

DECLARE_PROCESSOR(VertexAnaProcessor);
