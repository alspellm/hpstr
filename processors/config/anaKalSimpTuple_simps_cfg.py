import HpstrConf
import sys
import os
import baseConfig as base

base.parser.add_argument("-w", "--tracking", type=str, dest="tracking",
                  help="Which tracking to use to make plots", metavar="tracking", default="BOTH")
base.parser.add_argument("-V", "--splitVolume", type=int, dest="splitVolume",
                  help="Require positron in Top and Bottom", metavar="splitVolume", default=0)
base.parser.add_argument("-R", "--region", type=str, dest="region",
                  help="Signal Region (SR) or Control Region (CR)", metavar="region", default="BOTH")
options = base.parser.parse_args()

# Use the input file to set the output file name
infile = options.inFilename
outfile = options.outFilename

#outfile = outfile.split(".root")[0]+"_"+options.tracking+".root"
outfile = outfile.split(".root")[0]+".root"

print('Input file: %s' % infile)
print('Output file: %s' % outfile)

p = HpstrConf.Process()

p.run_mode = 1
#p.max_events = 1000

# Library containing processors
p.add_library("libprocessors")

###############################
#          Processors         #
###############################

recoana_kf = HpstrConf.Processor('vtxana_kf', 'VertexAnaProcessor')
recoana_gbl = HpstrConf.Processor('vtxana_gbl', 'VertexAnaProcessor')
mcana =  HpstrConf.Processor('mcpartana', 'MCAnaProcessor')
###############################
#   Processor Configuration   #
###############################
#RecoHitAna
recoana_kf.parameters["anaName"] = "vtxana_kf"
recoana_kf.parameters["trkColl"] = "KalmanFullTracks"
recoana_kf.parameters["tsColl"] = "TSData"
recoana_kf.parameters["vtxColl"] = "UnconstrainedV0Vertices_KF"
recoana_kf.parameters["mcColl"]  = "MCParticle"
recoana_kf.parameters["hitColl"] = "SiClusters"
recoana_kf.parameters["ecalColl"] = "RecoEcalClusters_KF"
recoana_kf.parameters["histoCfg"] = os.environ['HPSTR_BASE']+"/analysis/plotconfigs/tracking/vtxAnalysis.json"
recoana_kf.parameters["mcHistoCfg"] = os.environ['HPSTR_BASE']+'/analysis/plotconfigs/mc/basicMC.json'
#####
recoana_kf.parameters["beamE"] = base.beamE[str(options.year)]
recoana_kf.parameters["isData"] = options.isData
recoana_kf.parameters["analysis"] = options.analysis
recoana_kf.parameters["debug"] = 0
CalTimeOffset=-999

if (options.isData==1):
    CalTimeOffset=55.
    print("Running on data file: Setting CalTimeOffset %d"  % CalTimeOffset)

elif (options.isData==0):
    CalTimeOffset=43.
    print("Running on MC file: Setting CalTimeOffset %d"  % CalTimeOffset)
else:
    print("Specify which type of ntuple you are running on: -t 1 [for Data] / -t 0 [for MC]")

recoana_kf.parameters["CalTimeOffset"]=CalTimeOffset
#Region definitions

#RegionPath=os.environ['HPSTR_BASE']+"/analysis/selections/simps/nsigma_matcher/"
RegionPath=os.environ['HPSTR_BASE']+"/analysis/selections/simps/mindist_matcher/"

recoana_kf.parameters["regionDefinitions"] = [RegionPath+'vertexSelection_MOUSEcuts_only.json',
                                              RegionPath+'vertexSelection_bumphunt.json',
                                              RegionPath+'vertexSelection_loose.json',
                                              RegionPath+'vertexSelection_bumphunt_Tight.json',
                                              RegionPath+'vertexSelection_bumphunt_Tight_nvtx.json',
                                              RegionPath+'vertexSelection_bumphunt_Tight_L1.json',
                                              RegionPath+'vertexSelection_bumphunt_L1.json',
                                              RegionPath+'vertexSelection_bumphunt_MOUSEcuts.json']
#RecoHitAna
recoana_gbl.parameters = recoana_kf.parameters.copy()
recoana_gbl.parameters["anaName"] = "vtxana_gbl"
recoana_gbl.parameters["vtxColl"] = "UnconstrainedV0Vertices"
recoana_gbl.parameters["tsColl"]   = "TSData"
recoana_gbl.parameters["hitColl"] = "RotatedHelicalTrackHits"
recoana_gbl.parameters["trkColl"] = "GBLTracks"
recoana_gbl.parameters["mcColl"]  = "MCParticle"
recoana_gbl.parameters["ecalColl"] = "RecoEcalClusters_GBL" #NEW
#recoana_gbl.parameters["ecalColl"] = "RecoEcalClusters" #OLD
recoana_gbl.parameters["histoCfg"] = os.environ['HPSTR_BASE']+"/analysis/plotconfigs/tracking/vtxAnalysis.json"
recoana_gbl.parameters["mcHistoCfg"] = os.environ['HPSTR_BASE']+'/analysis/plotconfigs/mc/basicMC.json'
#####
recoana_gbl.parameters["beamE"] = base.beamE[str(options.year)]
recoana_gbl.parameters["isData"] = options.isData
recoana_gbl.parameters["analysis"] = options.analysis
recoana_gbl.parameters["debug"] = 0
recoana_gbl.parameters["CalTimeOffset"]=CalTimeOffset

#RegionPath=os.environ['HPSTR_BASE']+"/analysis/selections/simps/nsigma_matcher/"
RegionPath=os.environ['HPSTR_BASE']+"/analysis/selections/simps/mindist_matcher/"

recoana_gbl.parameters["regionDefinitions"] = [RegionPath+'vertexSelection_MOUSEcuts_only.json',
                                              RegionPath+'vertexSelection_bumphunt.json',
                                              RegionPath+'vertexSelection_loose.json',
                                              RegionPath+'vertexSelection_bumphunt_Tight.json',
                                              RegionPath+'vertexSelection_bumphunt_Tight_nvtx.json',
                                              RegionPath+'vertexSelection_bumphunt_Tight_L1.json',
                                              RegionPath+'vertexSelection_bumphunt_L1.json',
                                              RegionPath+'vertexSelection_bumphunt_MOUSEcuts.json']


#MCParticleAna
mcana.parameters["debug"] = 0
mcana.parameters["anaName"] = "mcAna"
mcana.parameters["partColl"] = "MCParticle"
mcana.parameters["trkrHitColl"] = "TrackerHits"
mcana.parameters["ecalHitColl"] = "EcalHits"
mcana.parameters["histCfg"] = os.environ['HPSTR_BASE']+'/analysis/plotconfigs/mc/basicMC.json'

#
#    RegionPath+'ESumCR.json',
#    RegionPath+'TightNoSharedL0.json',
#    RegionPath+'TightNoShared.json',

# Sequence which the processors will run.
#p.sequence = [recoana_kf,recoana_gbl]
if (options.tracking == "KF"):
    print("Run KalmanFullTracks analysis")
    p.sequence = [recoana_kf] #,mcana]
elif (options.tracking == "GBL"):
    print("Run GBL analysis")
    p.sequence = [recoana_gbl]#,mcana]
elif (options.tracking == "BOTH"):
    print("Run GBL and KF analysis")
    p.sequence = [recoana_gbl, recoana_kf]#,mcana]
else :
    print ("ERROR::Need to specify which tracks KF or GBL")
    exit(1)

if (options.nevents > 0):
    p.max_events = options.nevents

p.input_files=infile
p.output_files = [outfile]

p.printProcess()
