#include "BlFitHistos.h"

BlFitHistos::BlFitHistos(const std::string& inputName) {
    m_name = inputName;
    mmapper_ = new ModuleMapper();
}

BlFitHistos::~BlFitHistos() {
    for (std::map<std::string, TGraphErrors*>::iterator it = baselineGraphs.begin();
                it!=baselineGraphs.end(); ++it) {
            if (it->second) {
                delete (it->second);
                it->second = nullptr;
            }
        }
        baselineGraphs.clear();
}



void BlFitHistos::Chi2GausFit( HistoManager* inputHistos_, int nPointsDer_,int rebin_,int xmin_, int minStats_, FlatTupleMaker* flat_tuple_) {
     

    //Loop over all 2D histogram keys found in input ROOT File
    for (std::vector<std::string>::iterator jj = inputHistos_->histos2dNamesfromTFile.begin();
            jj != inputHistos_->histos2dNamesfromTFile.end(); ++jj)
    {
        //Get 2D Histograms from input ROOT File using keys, and rebin based on configuration setting "rebin"
        TH2F* histo_hh = inputHistos_->get2dHisto(*jj);
        histo_hh->RebinY(rebin_);
        histo_hh->Write();
        //Extract histogram key by removing _hh from key
        std::string SvtAna2DHisto_key = histo_hh->GetName();
        SvtAna2DHisto_key.erase(SvtAna2DHisto_key.end()-3,SvtAna2DHisto_key.end());

        //Perform fitting procedure over all channels on a sensor
        for(int cc=0; cc < 640 ; ++cc) 
        {
            double TFRE = 0.0;
            //Set Channel and Hybrid information and paramaters in the flat tuple
            flat_tuple_->setVariableValue("SvtAna2DHisto_key", SvtAna2DHisto_key);
            flat_tuple_->setVariableValue("channel", cc);
            flat_tuple_->setVariableValue("minbinThresh",(double)xmin_);
            flat_tuple_->setVariableValue("minStats", (double)minStats_);
            flat_tuple_->setVariableValue("rebin", (double)rebin_);


            std::cout <<"Hybrid: "<< SvtAna2DHisto_key << " channel #" << cc << std::endl;

            //Get YProjection (1D Channel Histo) from 2D Histogram 
            TH1D* projy_h = histo_hh->ProjectionY(Form("%s_proY_ch%i",SvtAna2DHisto_key.c_str(),cc),
                    cc+1,cc+1,"e");
            projy_h->SetTitle(Form("%s_proY_ch%i",SvtAna2DHisto_key.c_str(),cc));

            //Check number of entries and RMS of channel
            flat_tuple_->setVariableValue("n_entries", projy_h->GetEntries());
            double chRMS = projy_h->GetRMS();
            //If rms is > 400, set variable "noisy" to 1 for true
            flat_tuple_->setVariableValue("rms", chRMS);
            if(chRMS > 400) 
            {
                flat_tuple_->setVariableValue("noisy", 1.0);
            }
            else 
            {
                flat_tuple_->setVariableValue("noisy", 0.0);
            }


            //The fit window minimum range is found to be the first bin containing a number of
            //entries >= some fraction of the maximum number of entries in the 1D histogram.
            //If that fraction is less than the configurable minimum bin threshold value, instead
            //find the first bin above the minimum bin threshold. If no bin is found, skip fit
            double maxbin = projy_h->GetBinContent(projy_h->GetMaximumBin());
            int firstbin = 0;
            double frac = 0.15;
            std::cout << "[BlFitHistos] Find first bin above: " << frac*maxbin << std::endl;
            if (frac*maxbin <= xmin_)
            {
                firstbin = projy_h->FindFirstBinAbove((double)xmin_,1);
                std::cout << "[BlFitHistos] First bin: " << firstbin << std::endl;
            }
            else 
            {
                firstbin = projy_h->FindFirstBinAbove((double)frac*maxbin,1);
                std::cout << "[BlFitHistos] First bin: " << firstbin << std::endl;
            }
            
            //If channel does not have the minimum statistics required, set all variables to -9999.9
            //and skip the fit procedure on this channel
            std::cout << " skip channel?" << std::endl;
            if (firstbin == -1 or projy_h->GetEntries() < minStats_ ) 
            {
                flat_tuple_->setVariableValue("BlFitMean", -9999.9);
                flat_tuple_->setVariableValue("BlFitSigma", -9999.9);
                flat_tuple_->setVariableValue("BlFitNorm", -9999.9);
                flat_tuple_->setVariableValue("BlFitRangeLower", -9999.9);
                flat_tuple_->setVariableValue("BlFitRangeUpper", -9999.9);
                flat_tuple_->setVariableValue("BlFitChi2", -9999.9);
                flat_tuple_->setVariableValue("BlFitNdf", -9999.9);

                flat_tuple_->setVariableValue("ogxmax", -9999.9);
                flat_tuple_->setVariableValue("ogxmin", -9999.9);

                flat_tuple_->addToVector("iterMean", -9999.9);
                flat_tuple_->addToVector("iterChi2NDF_2der",-9999.9);
                flat_tuple_->addToVector("iterChi2NDF_1der",-9999.9);
                flat_tuple_->addToVector("iterChi2NDF_derRange",-9999.9);
                flat_tuple_->addToVector("iterChi2NDF",-9999.9);
                flat_tuple_->addToVector("ratio2derChi2", -9999.9);
                flat_tuple_->addToVector("iterFitRangeEnd", -9999.9);

                flat_tuple_->setVariableValue("lowdaq", -9999.9);
                flat_tuple_->setVariableValue("TFitResultError", -9999.9);

                flat_tuple_->setVariableValue("minbinFail",1.0);
                std::cout << "Not enough stats in channel " << cc << std::endl;
                flat_tuple_->fill();
                continue;
            }
            flat_tuple_->setVariableValue("minbinFail",0.0);

            //xmin is the start of the fit window. iterxmax will initially be iteratively fit and
            //then increased until some maximum allowed value, or until the chi2/Ndf > 100
            int iter = 0;
            double xmin = projy_h->GetBinLowEdge(firstbin);
            double binwidth = projy_h->GetBinWidth(firstbin);
            double iterxmax = xmin + 20.0*binwidth;
            std::cout << "[BlFitHistos] Initial xmin: " << xmin << std::endl;
            std::cout << "[BlFitHistos] Initial iterxmax: " << iterxmax << std::endl;
            
            //Define all iterative and final fit parameters
            std::vector<double> amp,mean,sigma,chi2,const_err,sigma_err,
                fit_range_end,chi2_NDF,chi2_2D,chi2_1D,der2chi2,der2chi2R;
            std::vector<int> NDF;


            //Baseline signals are composed of a gaussian baseline, followed by a landau pile-up
            //distribution. It is found that at the boundary of these two distributions, the
            //second derivative of the Chi2/Ndf of a fit is maximized. 
            //Therefore, in order to determine the correct fit window, i.e. xmin and xmax, that
            //will fit a gaussian over just the baseline distribution, the fit works as follows:
            //Iterate xmax across the ADC distribution. At each new value of iterxmax, perform
            //a fit on the histogram. Add Chi2/Ndf for the iterxmax of each fit to a vector.
            //Locate the iterxmax value that corresponds to the maximum Chi2/Ndf 2nd derivative. This
            //is where the fit window should end (xmax) 

            std::cout << "[BlFitHistos] Starting Iterative Fit Procedure" << std::endl;
            double currentChi2 = 0.0;
            while(iterxmax < 6800.0 && currentChi2 < 100.0 || iter < 10)
            {
                //std::cout << "[BlFitHistos] Iterating iterxmax: " << iterxmax << std::endl;
                TFitResultPtr cc_fit = projy_h->Fit("gaus", "QRES", "", xmin, iterxmax);
                if(cc_fit->Ndf() == 0)
                {
                    //std::cout << "[BlFitHistos] WARNING: Ndf = 0. SKIP CH" << firstbin << std::endl;
                    iterxmax = iterxmax + binwidth;
                    continue;
                }

                //Return iterative fit parameters and add to flat tuple
                flat_tuple_->addToVector("iterMean", cc_fit->GetParams()[1]);
                chi2.push_back(cc_fit->Chi2());
                NDF.push_back(cc_fit->Ndf());
                chi2_NDF.push_back(chi2.at(iter)/NDF.at(iter));
                flat_tuple_->addToVector("iterChi2NDF", chi2_NDF.at(iter));
                fit_range_end.push_back(iterxmax);
                flat_tuple_->addToVector("iterFitRangeEnd", iterxmax);

                //Increase fit range window and iterate the fit
                iterxmax = iterxmax + binwidth;
                iter++;
                currentChi2 = cc_fit->Chi2()/cc_fit->Ndf();
            }

            //Calculate the 2nd derivative of chi2 by taking slope on either side of one point
            //and taking the difference between the two. The xposition of the maximum Chi2 2nd 
            //derivative is used to define the maximum fit range x-position

            for(int i = nPointsDer_; i < chi2_NDF.size() - nPointsDer_; i++) 
            {
                double derForward=(chi2_NDF.at(i+nPointsDer_)-chi2_NDF.at(i))/(nPointsDer_*binwidth);
                double derBack = (chi2_NDF.at(i)-chi2_NDF.at(i-nPointsDer_))/(nPointsDer_*binwidth);
                double der = (derForward+derBack)/2.0;
                double der2 = (derForward - derBack); //(nPointsDer_*binwidth);
                if(der2 == der2) 
                {
                    chi2_2D.push_back(der2);
                }
                else
                {
                    chi2_2D.push_back(-9999.9);
                }

                if(der == der) 
                {
                    chi2_1D.push_back(der);
                }
                else
                {
                    chi2_1D.push_back(-9999.9);
                }
                flat_tuple_->addToVector("iterChi2NDF_2der",der2);
                flat_tuple_->addToVector("iterChi2NDF_1der",der);
            }

            //NO LONGER USED...
            //Take ratio of chi2/NDF 2nd derivative and chi2/NDF. ADC value of this maximum will be
            //the maximum fit range.
            std::vector<double> ratio2derChi2;
            for(int i=0; i < chi2_2D.size(); ++i)
            {
                ratio2derChi2.push_back(chi2_2D.at(i)/chi2_NDF.at(i+nPointsDer_));
                flat_tuple_->addToVector("ratio2derChi2",ratio2derChi2.at(i));
            }

            //FOLLOWING SECTION NEEDS TO BE UPDATED TO A CORRECT METHOD (advance iterator)
            //Create subrange for chi2_2D that accounts for derivative requiring <n> points prior  
            //and post point of interest
            double minimum_thresh = 1.;
            std::vector<double>::const_iterator first = fit_range_end.begin()+nPointsDer_;
            std::vector<double>::const_iterator last=fit_range_end.begin()+nPointsDer_+chi2_2D.size();
            std::vector<double> chi2_2D_range(first,last);
            
            //Find maximum Chi2/Ndf
            first = chi2_NDF.begin()+nPointsDer_;
            last = chi2_NDF.begin() + nPointsDer_ + chi2_2D.size();
            std::vector<double> chi2SubRange(first,last);
            double chi2max = *std::max_element(chi2SubRange.begin(), chi2SubRange.end());
            int chi2maxIndex = std::max_element(chi2SubRange.begin(), chi2SubRange.end()) - chi2SubRange.begin();

            //Find maximum chi2 second derivative
            double chi2_2D_max = *std::max_element(chi2_2D.begin(), chi2_2D.end());
            int chi2_2D_maxIndex = std::max_element(chi2_2D.begin(), chi2_2D.end()) - chi2_2D.begin();
            double chi2_2D_xmax = chi2_2D_range.at(chi2_2D_maxIndex);

            //xmax defines the end of the fit window over the gaussian distribution. Ideally this
            //position occurs at the boundary of the baseline and landau pile-up distributions
            double xmax = chi2_2D_xmax;
            std::cout << "[BlFitHistos] First xmax:" << xmax << std::endl;

            //At this stage, there are a few possible outcomes from the current fit:
            //1. xmin and xmax correctly correspond to the range of the baseline distribution.
            //The fit is a nice gaussian that closely follows the data to some degree.
            //2. Chosen xmax does not correlate to the end of the baseline, and extends 
            //into the pile-up signal.The fit spans across the entire distribution and does not
            //at all resemble the baseline.
            //3. Chosen xmax does not correlate to to to the end of the basline, and occurs too early
            //in the distribution. The fit only captures some small fraction of the baseline and 
            //returns a gaussian with a large width.
            //4. LowDaq-Threshold cuts off the baseline, so the fit returns a gaussian with 
            //a large width, but that accurately represents the baseline signal. These are flagged
            //as "lowdaq" channels.
            //If fit mean + N*Sigma > xmax, channel has low daq threshold *or* bad fit
   
            
            //When using maximum Chi2/Ndf 2nd derivative to locate the fit window range, there
            //are sometimes anomalous large spikes in these vales  with no correlation to the 
            //expected fit window end. If a spike exists, the position that it occurs will be
            //selected as xmax, the end of the fit window range, causing a bad fit.
            //It is found that these spikes often occur after xmax values beyond the maximum Chi2/Ndf
            //The position of maximum Chi2/Ndf is then used as a cut. Only xmax values that occur
            //before maximum Chi2/Ndf of the iterative fits are kept.
            
            std::vector<double>::const_iterator begin = chi2_2D.begin();
            std::vector<double>::const_iterator end  = chi2_2D.begin()+(chi2maxIndex);
            std::vector<double> cutrange(begin,end);
            double cut_maxIndex = std::max_element(cutrange.begin(),cutrange.end())-cutrange.begin();

            double cutxmax = chi2_2D_range.at(cut_maxIndex);
            double cutxmax_1 = cutxmax;
            std::cout << "[BlFitHistos] xmax after Chi2Max cut: " << cutxmax << std::endl;


            //There still exist cases where the anomalous Chi2/Ndf 2nd derivative spikes occur
            //before the maximum Chi2/Ndf position.
            //To remove these cases as well, the fit is performed using the xmax value after
            //the Chi2/Ndf max cut is applied. There are a few cases at play here:
            //Case 1: The Chi2 max cut had an effect, it changed the location of xmax, meaning
            //it successfully removed an anamolous spike
            //Case 2: The Chi2 max cut had no effect. 
            //Case 2.a: Cut had no effect, and there is no spike. 
            //Case 2.b: Cut had no effect, and spike is not removed. This needs addressing.
     
            //Perform fit using the current location of xmax
            double chi2_1;
            double sigma_1;
            TFitResultPtr fit = projy_h->Fit("gaus", "QRES", "", xmin, cutxmax_1);
            if ( fit == -1) 
            {
                TFRE = TFRE + 1.0;
                flat_tuple_->setVariableValue("TFitResultError",TFRE);
                continue;
            }
            else
            {
                chi2_1 = fit->Chi2()/fit->Ndf();
                sigma_1 = fit->GetParams()[2];
            }
            
            //If Case 1, keep current xmax position.
            //If Case 2: The Chi2 max cut had no effect on which iterxmax position is selected,
            //meaning that there is either no spike to be removed, or that the Chi2 max cut failed 
            //to remove the spike.
            //To ensure that any spike is removed,the iterxmax cut location is iterated backwards
            //until it is found that the position of maximum chi2/Ndf 2nd derivative changes.
            //The histogram is then fit using window sizes defined by the original xmax position,
            //and the newly decreased xmax position. If the Chi2/Ndf of the fit does not improve
            //by some defined percentage, keep the original xmax value. If it does improve by
            //some significant amount, this is indicative of the new xmax position providing
            //a better fit.

            double cutxmax_2;
            double chi2_2;
            double sigma_2;
            bool cut2=false;

            if(cutxmax >= chi2_2D_xmax) 
            {
                cut2 = true;
                while (cutxmax >= chi2_2D_xmax && chi2maxIndex > 0)
                {
                    begin = chi2_2D.begin();
                    end  = chi2_2D.begin()+(chi2maxIndex);
                    std::vector<double> cutrange(begin,end);
                    cut_maxIndex = std::max_element(cutrange.begin(),cutrange.end())-cutrange.begin();
                    cutxmax = chi2_2D_range.at(cut_maxIndex);
                    cutxmax_2 = cutxmax;
                    chi2maxIndex = chi2maxIndex - 1;
                }
                fit = projy_h->Fit("gaus", "QRES", "", xmin, cutxmax_2);
                if ( fit == -1) 
                {
                    TFRE = TFRE + 1.0;
                    flat_tuple_->setVariableValue("TFitResultError",TFRE);
                    continue;
                }
                else
                {
                    chi2_2 = fit->Chi2()/fit->Ndf();
                    sigma_2 = fit->GetParams()[2];
                }
            }
            std::cout << "[BlFitHistos] chi2_1: " << chi2_1 << std::endl;
            std::cout << "[BlFitHistos] chi2_2:  " << chi2_2 << std::endl;
            std::cout << "[BlFitHistos] sigma_1: " << sigma_1 << std::endl;
            std::cout << "[BlFitHistos] sigma_2:  " << sigma_2 << std::endl;
            
            //If the second cut was applied, keep the second fit if the chi2/ndf improves 
            //significantly, but also only if the second fit did not return a partial gaussian
            //(meaning half or less of the full gaussian shape)
            double N = 0.2;
            if (cut2 == true) 
            {
                if (chi2_2 < 0.8*chi2_1 &&  cutxmax_2 > fit->GetParams()[1] + N*fit->GetParams()[2])
                {
                    xmax = cutxmax_2;
                }
                else 
                {
                    xmax = cutxmax_1;
                }
            }
            else 
            {
                xmax = cutxmax_1;
            }
        
            std::cout << "[BlFitHistos] xmax after spike removal: " << xmax << std::endl;
            //Variables for examining fit window distributions in event of issues
            flat_tuple_->setVariableValue("ogxmax", xmax);
            flat_tuple_->setVariableValue("ogxmin", xmin);

            //Refit the data using the new fit window range
            const double* fitparams;
            fit = projy_h->Fit("gaus", "QRES", "", xmin, xmax);
            if ( fit == -1) 
            {
                TFRE = TFRE + 1.0;
                flat_tuple_->setVariableValue("TFitResultError",TFRE);
                continue;
            }
            else
            {
                fitparams = fit->GetParams();
            }

            //LOWDAQ THRESHOLD FLAG
            //If fit mean + N*Sigma > xmax, channel has low daq threshold *or* bad fit
            //

            double maxvalavg = 0;
            for( int n=1; n < 11; n++) 
            {
               maxvalavg = maxvalavg + projy_h->GetBinContent(projy_h->GetMaximumBin()+n);
               if ( n == 10)
               {
                    maxvalavg = maxvalavg/10.0;
               }
            }

            //If a channel matches any of the below criteria:
            //1. If xmax or xmin occur within N*Sigma of the fit, indicative of a partial fit or a
            //lowDaq-threshold channel
            //2. If the maximum bin in a distribution occurs outside of the fit by some level of 
            //sigma, AND the average around the maximum bin is greater than the mean of the fit, the
            //channel is lowDaq but fitting an improper gaussian over a messy distribution. 
            //(VISUAL EXAMPLES OF DISTRIBUTIONS NEEDED TO UNDERSTAND)
                        
            if ((fitparams[1] + N*fitparams[2] > xmax) || (fitparams[1] - N*fitparams[2] < xmin) || ( (projy_h->GetBinLowEdge(projy_h->GetMaximumBin()) > fitparams[1] + 0.1*fitparams[2]) && (maxvalavg > fitparams[0]))) 
            {
                //Fix xmax to the location of the maximum bin value, which should generally
                //be the required end of the fit window for the baseline distribution.
                //Iterate xmin up until the fit conforms to the baseline.
                double tempxmin = xmin;
                xmax = projy_h->GetBinLowEdge(projy_h->GetMaximumBin()) - 2*binwidth;
                double ogChi2Ndf = fit->Chi2()/fit->Ndf();
                if (xmax-xmin < 300)
                {
                    tempxmin = xmax - 300;
                    
                }
                //Minimum Fit window size seems to be necessary to avoid empty TFitResultPtr...?
                //200 chosen based on distribution of window sizes, but may need adjusting...
                while(xmin < (xmax - 200.0))
                {
                    std::cout << "[BlFitHistos] lowdaq xmin: " << tempxmin << std::endl;
                    tempxmin = tempxmin + binwidth;
                    TFitResultPtr fit = projy_h->Fit("gaus", "QRES", "", tempxmin, xmax);
                    if ( fit == -1) 
                    {
                        TFRE = TFRE + 1.0;
                        flat_tuple_->setVariableValue("TFitResultError",TFRE);
                        continue;
                    }
                    else
                    {
                        if ( fit->Chi2()/fit->Ndf() < 0.8*ogChi2Ndf )
                        {
                            ogChi2Ndf = fit->Chi2()/fit->Ndf();
                            xmin = tempxmin;
                            std::cout << "Increasing xmin to: " << xmax << std::endl;
                        }
                       else 
                       {
                            break;
                       }
                    }

                std::cout << "chi2/ndf is " << fit->Chi2()/fit->Ndf() << std::endl;

                }
            }


            //After attempting to refit baseline above, if lowDaq criteria is still met, label
            //channel as lowDaq and stop trying to refit.
            if ( fitparams[1] + N*fitparams[2] > xmax || fitparams[1] - N*fitparams[2] < xmin  )
            {
                flat_tuple_->setVariableValue("lowdaq", 1.0);
            }
            else 
            {
                double ogChi2Ndf = fit->Chi2()/fit->Ndf();
                std::cout << "original Chi2/Ndf is: " << ogChi2Ndf << std::endl;
                while ( xmax > (xmin + 100))
                {
                   double tempxmax = xmax;
                   tempxmax = tempxmax - 5*binwidth; 
                   TFitResultPtr newfit = projy_h->Fit("gaus", "QRES", "", xmin, tempxmax);
                    if ( newfit == -1) 
                    {
                        TFRE = TFRE + 1.0;
                        flat_tuple_->setVariableValue("TFitResultError",TFRE);
                        continue;
                    }
                   else 
                   {
                       std::cout << "New Chi2/Ndf is: " << newfit->Chi2()/newfit->Ndf() << std::endl;
                       if ( newfit->Chi2()/newfit->Ndf() < 0.8*ogChi2Ndf )
                       {
                            ogChi2Ndf = newfit->Chi2()/newfit->Ndf();
                            xmax = tempxmax;
                            std::cout << "Stepping back xmax to: " << xmax << std::endl;
                            std::cout << "xmin is: " << xmin << std::endl;
                       }
                       else 
                       {
                            break;
                       }
                   }
                }

                flat_tuple_->setVariableValue("lowdaq", 0.0);
            }

            fit = projy_h->Fit("gaus", "QRES", "", xmin, xmax);
            if ( fit == -1) 
            {
                TFRE = TFRE + 1.0;
                flat_tuple_->setVariableValue("TFitResultError",TFRE);
                continue;
            }
            else
            {
                flat_tuple_->setVariableValue("TFitResultError",TFRE);
                fitparams = fit->GetParams();
            }

            //Store Final Fit Parameters in the flat tuple
            //std::cout << "New Xmax-Xmin = " << xmax - xmin << std::endl;
            flat_tuple_->setVariableValue("BlFitMean", fitparams[1]);
            flat_tuple_->setVariableValue("BlFitSigma", fitparams[2]);
            flat_tuple_->setVariableValue("BlFitNorm", fitparams[0]);
            flat_tuple_->setVariableValue("BlFitRangeLower", (double)xmin);
            flat_tuple_->setVariableValue("BlFitRangeUpper", (double)xmax);
            flat_tuple_->setVariableValue("BlFitChi2", fit->Chi2());
            flat_tuple_->setVariableValue("BlFitNdf", fit->Ndf());
            //std::cout << "Final Chi2/Ndf is " << fit->Chi2()/fit->Ndf();
            for(int i=0; i < chi2_2D_range.size(); ++i) 
            {
                flat_tuple_->addToVector("iterChi2NDF_derRange",chi2_2D_range.at(i));
            }

        flat_tuple_->fill();
        delete projy_h;
        }
    }

}

