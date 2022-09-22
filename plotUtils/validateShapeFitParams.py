#!/usr/bin/python3
import ROOT as r
import numpy as np
import argparse

parser=argparse.ArgumentParser(description="")
parser.add_argument("-i", type=str, dest="inFilename", help="Input File",default="")
parser.add_argument("-o", type=str, dest="outFilename", help="output root file",default="pulseFit_ana.root")
options = parser.parse_args()

def getFitTuple(inFile):

    hwTag, channel, svtid, layer, module, t0, tau1, tau2, amp, t0err, tau1err, tau2err, amperr, chi2, ndf = ([] for i in range(15))
    tree = inFile.fits
    for e in tree:
        print(e.svtid)


###############################################################################################

inFile = r.TFile(options.inFilename,"READ")
outFile = r.TFile(options.outFilename,"RECREATE")
tree = inFile.fits

####### Histograms

slim_tau12_hh = r.TH2F("slim_sensor_tau1_v_2","slim_sensor_tau1_v_2;tau1;tau2",1000,0,100,1000,0,100)
slim_t0amp_hh = r.TH2F("slim_sensor_t0_v_amp","slim_sensor_t0_v_amp;t0 (ns);amp (adc)",60,-60,0,4000,0,400000)
slim_tau1_h = r.TH1F("slim_sensor_tau1","slim_sensor_tau1;svtid;tau1",5000,0,5000)
slim_tau2_h = r.TH1F("slim_sensor_tau2","slim_sensor_tau2;svtid;tau2",5000,0,5000)
slim_t0_h = r.TH1F("slim_sensor_t0","slim_sensor_t0;svtid;t0 (ns);",5000,0,5000)

thick_tau12_hh = r.TH2F("thick_sensor_tau1_v_2","thick_sensor_tau1_v_2;tau1;tau2",1000,0,100,1000,0,100)
thick_t0amp_hh = r.TH2F("thick_sensor_t0_v_amp","thick_sensor_t0_v_amp;t0 (ns);amp (adc)",60,-60,0,4000,0,400000)
thick_tau1_h = r.TH1F("thick_sensor_tau1","thick_sensor_tau1;svtid;tau1",25000,0,25000)
thick_tau2_h = r.TH1F("thick_sensor_tau2","thick_sensor_tau2;svtid;tau2",25000,0,25000)
thick_t0_h = r.TH1F("thick_sensor_t0","thick_sensor_t0;svtid;t0 (ns);",25000,0,25000)

all_tau12_hh = r.TH2F("all_sensor_tau1_v_2","all_sensor_tau1_v_2;tau1;tau2",1000,0,100,1000,0,100)
all_t0amp_hh = r.TH2F("all_sensor_t0_v_amp","all_sensor_t0_v_amp;t0 (ns);amp (adc)",60,-60,0,4000,0,400000)
all_tau1_h = r.TH1F("all_sensor_tau1","all_sensor_tau1;svtid;tau1",25000,0,25000)
all_tau2_h = r.TH1F("all_sensor_tau2","all_sensor_tau2;svtid;tau2",25000,0,25000)
all_t0_h = r.TH1F("all_sensor_t0","all_sensor_t0;svtid;t0 (ns);",25000,0,25000)

tau1err_hh = r.TH2F("all_sensor_tau1_v_err","all_sensor_tau1_v_err;tau1;err",100,0,100,100,0,100)
tau2err_hh = r.TH2F("all_sensor_tau2_v_err","all_sensor_tau2_v_err;tau2;err",100,0,100,100,0,100)
t0err_hh = r.TH2F("all_sensor_t0_v_err","all_sensor_t0_v_err;t0 (ns);err",60,-60,0,100,0,100)
amperr_hh = r.TH2F("all_sensor_amp_v_err","all_sensor_amp_v_err;amp (adc);err",4000,0,400000,1000,0,10000)

tau12err_hh = r.TH2F("all_sensor_tau1_v_2_err","all_sensor_tau1_v_2_err;tau1err;tau2err",100,0,100,100,0,100)
#Histos post error cuts

slim_tau12cut_hh = r.TH2F("slim_sensor_err_cut_tau1_v_2","slim_sensor_err_cut_tau1_v_2;tau1;tau2",1000,0,100,1000,0,100)
slim_t0ampcut_hh = r.TH2F("slim_sensor_err_cut_t0_v_amp","slim_sensor_err_cut_t0_v_amp;t0 (ns);amp (adc)",60,-60,0,4000,0,400000)
slim_tau1cut_h = r.TH1F("slim_sensor_err_cut_tau1","slim_sensor_err_cut_tau1;svtid;tau1",5000,0,5000)
slim_tau2cut_h = r.TH1F("slim_sensor_err_cut_tau2","slim_sensor_err_cut_tau2;svtid;tau2",5000,0,5000)
slim_t0cut_h = r.TH1F("slim_sensor_err_cut_t0","slim_sensor_err_cut_t0;svtid;t0 (ns);",5000,0,5000)

thick_tau12cut_hh = r.TH2F("thick_sensor_err_cut_tau1_v_2","thick_sensor_err_cut_tau1_v_2;tau1;tau2",1000,0,100,1000,0,100)
thick_t0ampcut_hh = r.TH2F("thick_sensor_err_cut_t0_v_amp","thick_sensor_err_cut_t0_v_amp;t0 (ns);amp (adc)",60,-60,0,4000,0,400000)
thick_tau1cut_h = r.TH1F("thick_sensor_err_cut_tau1","thick_sensor_err_cut_tau1;svtid;tau1",25000,0,25000)
thick_tau2cut_h = r.TH1F("thick_sensor_err_cut_tau2","thick_sensor_err_cut_tau2;svtid;tau2",25000,0,25000)
thick_t0cut_h = r.TH1F("thick_sensor_err_cut_t0","thick_sensor_err_cut_t0;svtid;t0 (ns);",25000,0,25000)

all_tau12cut_hh = r.TH2F("all_sensor_err_cut_tau1_v_2","all_sensor_err_cut_tau1_v_2;tau1;tau2",1000,0,100,1000,0,100)
all_t0ampcut_hh = r.TH2F("all_sensor_err_cut_t0_v_amp","all_sensor_err_cut_t0_v_amp;t0 (ns);amp (adc)",60,-60,0,4000,0,400000)
all_tau1cut_h = r.TH1F("all_sensor_err_cut_tau1","all_sensor_err_cut_tau1;svtid;tau1",25000,0,25000)
all_tau2cut_h = r.TH1F("all_sensor_err_cut_tau2","all_sensor_err_cut_tau2;svtid;tau2",25000,0,25000)
all_t0cut_h = r.TH1F("all_sensor_err_cut_t0","all_sensor_err_cut_t0;svtid;t0 (ns);",25000,0,25000)

tau1errcut_hh = r.TH2F("all_sensor_err_cut_tau1_v_err","all_sensor_err_cut_tau1_v_err;tau1;err",100,0,100,100,0,100)
tau2errcut_hh = r.TH2F("all_sensor_err_cut_tau2_v_err","all_sensor_err_cut_tau2_v_err;tau2;err",100,0,100,100,0,100)
t0errcut_hh = r.TH2F("all_sensor_err_cut_t0_v_err","all_sensor_err_cut_t0_v_err;t0 (ns);err",60,-60,0,100,0,100)


tau1chi2_hh = r.TH2F("all_sensor_tau1_v_chi2","all_sensor_tau1_v_chi2;tau1;chi2",100,0,100,1000,0,10)
tau2chi2_hh = r.TH2F("all_sensor_tau2_v_chi2","all_sensor_tau2_v_chi2;tau1;chi2",100,0,100,1000,0,10)
t0chi2_hh = r.TH2F("all_sensor_t0_v_chi2","all_sensor_t0_v_chi2;t0;chi2",100,0,100,1000,0,10)
chi2_h = r.TH1F("all_sensor_chi2","all_sensor_chi2",1000,0,10)
chi2_svtid_h = r.TH1F("all_sensor_svtid_chi2","all_sensor_svtid_chi2;svtid;chi2",25000,0,25000)


t0tau2_hh = r.TH2F("all_sensor_t0_v_tau2","all_sensor_t0_v_tau2;t0;tau2",1000,0,100,1000,0,100)
slim_t0tau2_hh = r.TH2F("slim_sensor_t0_v_tau2","slim_sensor_t0_v_tau2;t0;tau2",1000,0,100,1000,0,100)
thick_t0tau2_hh = r.TH2F("thick_sensor_t0_v_tau2","thick_sensor_t0_v_tau2;t0;tau2",1000,0,100,1000,0,100)

#mod group plots
tau1mod8_thick_hh = r.TH2F("thick_sensor_tau1_mod_8","thick_sensor_tau1_mod_8;mod_group;tau1",8,0,8,100,0,100)
tau2mod8_thick_hh = r.TH2F("thick_sensor_tau2_mod_8","thick_sensor_tau2_mod_8;mod_group;tau2",8,0,8,100,0,100)
t0mod8_thick_hh = r.TH2F("thick_sensor_t0_mod_8","thick_sensor_t0_mod_8;mod_group;t0",8,0,8,100,0,100)
chi2mod8_thick_hh = r.TH2F("thick_sensor_chi2_mod_8","thick_sensor_chi2_mod_8;mod_group;t0",8,0,8,250,0,500)

tau1mod8_slim_hh = r.TH2F("slim_sensor_tau1_mod_8","slim_sensor_tau1_mod_8;mod_group;tau1",8,0,8,100,0,100)
tau2mod8_slim_hh = r.TH2F("slim_sensor_tau2_mod_8","slim_sensor_tau2_mod_8;mod_group;tau2",8,0,8,100,0,100)
t0mod8_slim_hh = r.TH2F("slim_sensor_t0_mod_8","slim_sensor_t0_mod_8;mod_group;t0",8,0,8,100,0,100)
chi2mod8_slim_hh = r.TH2F("slim_sensor_chi2_mod_8","slim_sensor_chi2_mod_8;mod_group;t0",8,0,8,250,0,500)

#apv channel number plots
tau2_chan_thick_hh = r.TH2F("thick_sensor_channel_tau2","thick_sensor_channel_tau2;apv_channel;tau2",128,-0.5,127.5,1000,0,100)
tau1_chan_thick_hh = r.TH2F("thick_sensor_channel_tau1","thick_sensor_channel_tau1;apv_channel;tau1",128,-0.5,127.5,1000,0,100)

#######


hwTag, channel, svtid, layer, module, t0, tau1, tau2, amp, t0err, tau1err, tau2err, amperr, chi2, ndf, badPulse, nanfit = ([] for i in range(17))
badpulses = []
apv127s = []
nanfits = []
allfits = {}
goodfits = {}
nfits = 0

tau2largespikes = []
tau2lowspikes = []
tau1largespikes = []
tau1lowspikes = []

cutcount = 0
for e in tree:
    svtid = e.svtid
    layer = e.layer
    module = e.module
    channel = e.channel
    hwTag = e.hwTag
    t0 = e.t0
    tau1 = e.tau1
    tau2 = e.tau2
    amp = e.amp
    t0err = e.t0err
    tau1err = e.tau1err
    tau2err = e.tau2err
    amperr = e.amperr
    chi2 = e.chi2
    ndf = e.ndf
    badPulse = e.badPulse
    nanfit = e.nanfit
    apv127 = False

    fitparams = [svtid, layer, module, hwTag, t0, tau1, tau2, amp, t0err, tau1err, tau2err, amperr, chi2, ndf, badPulse, nanfit]
    allfits[svtid] = fitparams


    if badPulse:
        badpulses.append(svtid)
        continue
    if nanfit:
        nanfits.append(svtid)
        continue

    if svtid > 4095 and (channel%128 == 127 or channel%128 == 0):
        cutcount = cutcount + 1
        continue

    if channel%8 == 7:
        cutcount = cutcount + 1
        continue

    goodfits[svtid] = fitparams

    #error histos
    tau1err_hh.Fill(tau1,tau1err)
    tau2err_hh.Fill(tau2,tau2err)
    t0err_hh.Fill(t0,t0err)
    tau12err_hh.Fill(tau1err,tau2err)
    amperr_hh.Fill(amp,amperr)

    #all histos
    all_tau12_hh.Fill(tau1,tau2)
    all_t0amp_hh.Fill(t0,amp)
    all_tau1_h.SetBinContent(int(svtid)+1,tau1)
    all_tau2_h.SetBinContent(int(svtid)+1,tau2)
    all_t0_h.SetBinContent(int(svtid)+1,t0)
    t0tau2_hh.Fill(t0,tau2)


    if ndf > 0:
        tau1chi2_hh.Fill(tau1,chi2/ndf)
        tau2chi2_hh.Fill(tau2,chi2/ndf)
        t0chi2_hh.Fill(t0,chi2/ndf)
        chi2_h.Fill(chi2/ndf)
        chi2_svtid_h.SetBinContent(int(svtid)+1, chi2/ndf)

    #slim sensors
    if svtid < 4096:
        slim_tau12_hh.Fill(tau1,tau2)
        slim_t0amp_hh.Fill(t0,amp)
        slim_tau1_h.SetBinContent(int(svtid)+1,tau1)
        slim_tau2_h.SetBinContent(int(svtid)+1,tau2)
        slim_t0_h.SetBinContent(int(svtid)+1,t0)
        slim_t0tau2_hh.Fill(t0,tau2)

        #mod plots
        tau1mod8_slim_hh.Fill(channel%8,tau1)
        tau2mod8_slim_hh.Fill(channel%8,tau2)
        t0mod8_slim_hh.Fill(channel%8,t0)
        if ndf > 0:
            chi2mod8_slim_hh.Fill(channel%8,chi2/ndf)

    #thick sensors
    if svtid > 4095:
        thick_tau12_hh.Fill(tau1,tau2)
        thick_t0amp_hh.Fill(t0,amp)
        thick_tau1_h.SetBinContent(int(svtid)+1,tau1)
        thick_tau2_h.SetBinContent(int(svtid)+1,tau2)
        thick_t0_h.SetBinContent(int(svtid)+1,t0)
        thick_t0tau2_hh.Fill(t0,tau2)
        if ndf > 0:
            chi2mod8_thick_hh.Fill(channel%8,chi2/ndf)

        #mod plots
        tau1mod8_thick_hh.Fill(channel%8,tau1)
        tau2mod8_thick_hh.Fill(channel%8,tau2)
        t0mod8_thick_hh.Fill(channel%8,t0)

        apv_channel = channel%128
        tau2_chan_thick_hh.Fill(apv_channel,tau2) 
        tau1_chan_thick_hh.Fill(apv_channel,tau1) 

print("Total Number of fits: ", len(allfits))
print("Total Number of good fits: ", len(goodfits))
print("Number of NAN fits: ",len(nanfits))
print("Number of Bad Pulses: ", len(badpulses))
print("Bad pulse SVTIDS: ", badpulses)
print("NAN fits: ",(nanfits))
print("apv ch 127s skipped: ",len(apv127s))
print("cut channels ",cutcount)


#write histos
outFile.cd()

slim_tau12_hh.Write() 
slim_t0amp_hh.Write() 
slim_tau1_h.Write()
slim_tau2_h.Write()
slim_t0_h.Write()

thick_tau12_hh.Write()
thick_t0amp_hh.Write()
thick_tau1_h.Write()
thick_tau2_h.Write()
thick_t0_h.Write()

all_tau12_hh.Write()
all_t0amp_hh.Write()
all_tau1_h.Write()
all_tau2_h.Write()
all_t0_h.Write()

tau1chi2_hh.Write()
tau2chi2_hh.Write()
t0chi2_hh.Write()
chi2_h.Write()
t0tau2_hh.Write()
slim_t0tau2_hh.Write()
thick_t0tau2_hh.Write()
chi2_svtid_h.Write()

#mod plots
tau1mod8_thick_hh.Write()
tau2mod8_thick_hh.Write()
t0mod8_thick_hh.Write()
tau1mod8_slim_hh.Write()
tau2mod8_slim_hh.Write()
t0mod8_slim_hh.Write()
chi2mod8_slim_hh.Write()
chi2mod8_thick_hh.Write()

tau1err_hh.Write()
tau2err_hh.Write()
t0err_hh.Write()
tau12err_hh.Write()
amperr_hh.Write()

#chan plots
tau2_chan_thick_hh.Write()
tau1_chan_thick_hh.Write()


