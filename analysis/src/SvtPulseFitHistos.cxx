#include "SvtPulseFitHistos.h"
#include <math.h>
#include "TCanvas.h"

SvtPulseFitHistos::SvtPulseFitHistos(const std::string& inputName, ModuleMapper* mmapper) {
    m_name = inputName;
    mmapper_ = mmapper;
    svtid_map_ = mmapper_->buildChannelSvtIDMap();
}

SvtPulseFitHistos::~SvtPulseFitHistos() {
}

double SvtPulseFitHistos::GetHitTime(int sample_number, int cdel) {
    double time = 3.125*(8-cdel) + 25.0*sample_number;
    return time;
}

void SvtPulseFitHistos::buildRawSvtHitsTuple(std::vector<RawSvtHit*> *rawSvtHits_, FlatTupleMaker* rawhits_tup_) {
    
    int nhits = rawSvtHits_->size();
    for (int i = 0; i < nhits; i++)
    {
        RawSvtHit* rawSvtHit = rawSvtHits_->at(i);
        auto module = std::to_string(rawSvtHit->getModule());
        auto layer = std::to_string(rawSvtHit->getLayer());
        std::string hwTag= mmapper_->getHwFromSw("ly"+layer+"_m"+module);
        float channel = rawSvtHit->getStrip();
        int svtid = mmapper_->getSvtIDFromHWChannel(channel, hwTag, svtid_map_); 
        int calgroup = (int)channel%8;
        //Hard coded for now. Need to find way to determine this from data
        int cdel = 1;

        rawhits_tup_->setVariableValue("hwTag", hwTag);
        rawhits_tup_->setVariableValue("layer", layer);
        rawhits_tup_->setVariableValue("module", module);
        rawhits_tup_->setVariableValue("channel", channel);
        rawhits_tup_->setVariableValue("svtid", svtid);
        rawhits_tup_->setVariableValue("calgroup", calgroup);
        rawhits_tup_->setVariableValue("cdel", cdel);

        for (int ss = 0; ss < 6; ss++)
        {
            float adc = rawSvtHit->getADCs()[ss];
            rawhits_tup_->addToVector("adcs",adc);
        }

        rawhits_tup_->fill();
    }

}

void SvtPulseFitHistos::FillHistogramsByHw(std::vector<RawSvtHit*> *rawSvtHits_,float weight) {

    int nhits = rawSvtHits_->size();
    std::vector<std::string> hybridStrings={};
    std::string histokey;
    if(Event_number%10000 == 0) std::cout << "Event: " << Event_number 
        << " Number of RawSvtHits: " << nhits << std::endl;



    //Populates histograms for each hybrid channel
    for (int i = 0; i < nhits; i++)
    {
        RawSvtHit* rawSvtHit = rawSvtHits_->at(i);
        auto mod = std::to_string(rawSvtHit->getModule());
        auto lay = std::to_string(rawSvtHit->getLayer());
        std::string hwTag= mmapper_->getHwFromSw("ly"+lay+"_m"+mod);
        float channel = rawSvtHit->getStrip();
        int svtid = mmapper_->getSvtIDFromHWChannel(channel, hwTag, svtid_map_); 
            

        
        for (int ss = 0; ss < 6; ss++)
        {
            //histokey = hwTag + "_SvtHybrids_s"+std::to_string(ss)+"_hh";
            //histokey = hwTag + "_SvtHybrids_s0_hh";
            Fill2DHisto(histokey, 
                    (float)rawSvtHit->getStrip(),
                    (float)rawSvtHit->getADCs()[ss], 
                    weight);
        }
    }

    Event_number++;
}      
