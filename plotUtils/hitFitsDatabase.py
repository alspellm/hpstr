#!/usr/bin/python3
import ROOT as r
import numpy as np
import argparse
import csv
from operator import itemgetter

parser=argparse.ArgumentParser(description="")
parser.add_argument("-i", type=str, dest="inFilename", help="Input File",default="")
parser.add_argument("-o", type=str, dest="outFilename", help="output root file",default="pulseFit_ana.root")
parser.add_argument("-dbo", type=str, dest="databaseOutFile", help="db format output file",default="svt_shape_fit_params_2021.txt")

options = parser.parse_args()

def getFitTuple(inFile):

    hwTag, channel, svtid, layer, module, t0, tau1, tau2, amp, t0err, tau1err, tau2err, amperr, chi2, ndf = ([] for i in range(15))
    tree = inFile.fits
    for e in tree:
        print(e.svtid)

def getAmplitudeIntegralNorm(t1, t2):
    peak_t = 3.0*((t1* (t2**3))**(0.25))
    A = (t1**2)/((t1-t2)**3)
    B = (t1-t2)/(t1*t2)
    peakAmp = A * (np.exp(-peak_t / t1) - np.exp(-peak_t/t2) * (1 + peak_t * B + 0.5 * peak_t * peak_t *B*B))
    return peakAmp


###############################################################################################

inFile = r.TFile(options.inFilename,"READ")
tree = inFile.fits

hwTag, channel, svtid, layer, module, t0, tau1, tau2, amp, t0err, tau1err, tau2err, amperr, chi2, ndf, badPulse, nanfit = ([] for i in range(17))

badpulses = []
nanfits = []
allfits = {}
goodfits = {}
nfits = 0

dbout = []
existingids = {}

for e in tree:
    dbtuple = ()
    useneighbors = False

    svtid = e.svtid
    existingids[svtid] = 1
    layer = e.layer
    module = e.module
    channel = e.channel
    hwTag = e.hwTag
    
    t0 = e.t0
    tau1 = e.tau1
    tau2 = e.tau2
    amp = e.amp
    amp = amp*getAmplitudeIntegralNorm(tau1,tau2)
    t0err = e.t0err
    tau1err = e.tau1err
    tau2err = e.tau2err
    amperr = e.amperr
    chi2 = e.chi2
    ndf = e.ndf
    badPulse = e.badPulse
    nanfit = e.nanfit

    if svtid > 4095 and channel%128 == 127:
        badPulse = True

    if badPulse == True:
        useneighbors = True
    if nanfit == True:
        useneighbors = True

    dbtuple = (int(svtid), np.round(amp,2), round(t0,2), round(tau1,2), round(tau2,2), useneighbors)
    dbout.append(dbtuple)

dbout = sorted(dbout, key=itemgetter(0))
dboutF = options.databaseOutFile
for svtid in range(24576):
    if svtid not in existingids.keys():
        dbtuple = [int(svtid),0.0,0.0,0.0,0.0,True]
        dbout.append(dbtuple)

dbout = sorted(dbout, key=itemgetter(0))

with open(dboutF, 'w') as f:
    writer = csv.writer(f, delimiter = ',')
    #writer.writerow(['svt_channel_id', 'amplitude','t0','tp','tp2'])
    for i,e in enumerate(dbout):
        row = []
        nleft = 0
        nright = 0
        if e[5] == True:
            print("Found bad pulse on svtid ",e[0])
            #iterate left to find nearest channel with good pulse
            for ii in range(i-1,-1,-1):
                #print("Checking left of channel: svtid ", dbout[ii][0]) 
                nleft = nleft + 1
                if dbout[ii][5] == True:
                    continue
                else:
                    print("Use values of svtid ",dbout[ii][0], " for svtid: ", e[0])
                    row = [e[0], dbout[ii][1],dbout[ii][2],dbout[ii][3],dbout[ii][4]]
                    break

            #iterate right to find nearest channel with good pulse
            for ii in range(i+1,len(dbout)):
                #print("Checking right of channel: svtid ", dbout[ii][0]) 
                nright = nright + 1
                if dbout[ii][5] == True:
                    continue
                else:
                    if nright < nleft:
                        #print("Found closer channel to right")
                        print("Use values of svtid ",dbout[ii][0], " for svtid: ", e[0])
                        row = [e[0], dbout[ii][1],dbout[ii][2],dbout[ii][3],dbout[ii][4]]
                        break
                    else:
                        break
        #If pulse is good, just use values
        else:
            row = e[0:5]
        writer.writerow(row)
