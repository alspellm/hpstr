#include "Svt2DBlHistos.h"
#include <math.h>
#include "TCanvas.h"

Svt2DBlHistos::Svt2DBlHistos(const std::string& inputName, const std::string histoConfigFile) : HistoManager{inputName} {
    m_name_ = inputName;
    HistoManager{inputName};
    mmapper_ = new ModuleMapper();
    h_configs_ = histoConfigFile;
}

Svt2DBlHistos::~Svt2DBlHistos() {

    for (std::map<std::string, TGraphErrors*>::iterator it = baselineGraphs.begin(); 
            it!=baselineGraphs.end(); ++it) {
        if (it->second) {
            delete (it->second);
            it->second = nullptr;
        }
    }
    baselineGraphs.clear();
}

void Svt2DBlHistos::defineHybridHistos(std::vector<std::string> hybridNames, int histoDimension){


    if(histoDimension == 1){
        std::map<std::string,TH1F*>::iterator it = histos1d.begin();
        while (it!= histos1d.end()){
            std::string basename = it->first;
            TH1F* h = it->second;
            for(std::vector<std::string>::iterator it2 = hybridNames.begin(); it2 != hybridNames.end(); ++it2){
                std::string hybridName = *it2;
                std::string histoFullName = hybridName + "_" + basename;
                h->SetNameTitle(histoFullName.c_str(),histoFullName.c_str());
                histos1d.insert({histoFullName, h});
            }
            histos1d.erase(it++);
        }
        
    }

    else if(histoDimension ==2){
        std::map<std::string,TH2F*>::iterator it = histos2d.begin();
        while (it!= histos2d.end()){
            std::string basename = it->first;
            TH2F* hh = it->second;
            for(std::vector<std::string>::iterator it2 = hybridNames.begin(); it2 != hybridNames.end(); ++it2){
                std::string hybridName = *it2;
                std::string histoFullName = hybridName + "_" + basename;
                hh->SetNameTitle(histoFullName.c_str(),histoFullName.c_str());
                histos2d.insert({histoFullName, hh});
            }
            histos2d.erase(it++);
        }
        
    }
    
    /*
    else if(histoDimension == 2){
        std::map<std::string,TH2F*>::iterator it = histos2d.begin();
        while (it!= histos2d.end())
        {
            std::string basename = it->first;
            std::string histoFullName = hybridName + "_" + basename;
            TH2F* hh = it->second;
            hh->SetNameTitle(histoFullName.c_str(),histoFullName.c_str());
            histos2d.erase(basename);
            histos2d.insert({histoFullName, hh});
            std::cout << basename << " converted to: " << histos2d[histoFullName]->GetName() << std::endl;

            it++;

        }
    }
    */
       
}

void Svt2DBlHistos::buildHistos(){
    std::vector<std::string> hybridNames;
    mmapper_->getStrings(hybridNames);
    loadHistoConfig(h_configs_);
    //Define the template histograms
    DefineHistos();

    defineHybridHistos(hybridNames,1);
    defineHybridHistos(hybridNames,2);
         
    std::map<std::string,TH1F*>::iterator it = histos1d.begin();
    while (it!= histos1d.end())
    {
        std::cout << "TH1F " << it->first << "defined" << std::endl;
        it++;
    }
    std::map<std::string,TH2F*>::iterator it2 = histos2d.begin();
    while (it2!= histos2d.end())
    {
        std::cout << "TH2F " << it2->first << "defined" << std::endl;
        it2++;
    }


}


void Svt2DBlHistos::get2DHistoOccupancy(std::vector<std::string> histos2dNames) {


}



void Svt2DBlHistos::FillHistograms(std::vector<RawSvtHit*> *rawSvtHits_,float weight) {

    /*

    int nhits = rawSvtHits_->size();
    std::vector<std::string> hybridStrings={};
    std::string histokey;
    if(Event_number%1000 == 0) std::cout << "Event: " << Event_number 
        << " Number of RawSvtHits: " << nhits << std::endl;

    //Count the total number of hits each hybrid records per event
    int svtHybMulti[4][15] = {0};//build matrix that holds the different layer and moduel combinations
    for (int i = 0; i < nhits; i++)
    {
        RawSvtHit* rawSvtHit = rawSvtHits_->at(i);
        int mod = rawSvtHit->getModule();
        int lay = rawSvtHit->getLayer();
        svtHybMulti[mod][lay]++;

    }
    for (int i =0; i < 4; i++)
    {
        for (int j = 1; j < 15; j++)
        {
            if (!(j<9 && i>1))
            {   
                std::string swTag = mmapper_->getStringFromSw("ly"+std::to_string(j)+"_m"+std::to_string(i));
                hybridStrings.push_back(swTag);
                Fill1DHisto("hitN_"+swTag+"_h", svtHybMulti[i][j],weight);
            }
        }
    }

    Fill1DHisto("svtHitN_h", nhits,weight);
    //End of counting block

    //Populates histograms for each hybrid
    for (int i = 0; i < nhits; i++)
    {
        RawSvtHit* rawSvtHit = rawSvtHits_->at(i);
        auto mod = std::to_string(rawSvtHit->getModule());
        auto lay = std::to_string(rawSvtHit->getLayer());
        std::string swTag= mmapper_->getStringFromSw("ly"+lay+"_m"+mod);
        
        //Manually select which baselines (0 - 6) are included. THIS MUST MATCH THE JSON FILE!

        int ss = 0;
        
            //if(debug_ > 0) std::cout << "Filling Histogram with RawSvtHit" << std::endl; 
            histokey = "baseline"+std::to_string(ss)+"_"+swTag+"_hh";
                        Fill2DHisto(histokey, 
                    (float)rawSvtHit->getStrip(),
                    (float)rawSvtHit->getADCs()[ss], 
                    weight);
            //if(debug_ > 0) std::cout << "Histogram Filled" << std::endl; 
        

        ss = 3;
        
            //if(debug_ > 0) std::cout << "Filling Histogram with RawSvtHit" << std::endl; 
            histokey = "baseline"+std::to_string(ss)+"_"+swTag+"_hh";
                        Fill2DHisto(histokey, 
                    (float)rawSvtHit->getStrip(),
                    (float)rawSvtHit->getADCs()[ss], 
                    weight);
            //if(debug_ > 0) std::cout << "Histogram Filled" << std::endl; 
        
    }


            Event_number++;
            */
}      
