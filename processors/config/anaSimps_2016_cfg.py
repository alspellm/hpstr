import HpstrConf
import sys
import os
import baseConfig as base

base.parser.add_argument("-f", "--makeFlatTuple", type=int, dest="makeFlatTuple",
                         help="Make True to make vertex ana flat tuple", metavar="makeFlatTuple", default=1)
base.parser.add_argument("-r", "--isRadPDG", type=int, dest="isRadPDG",
                         help="Set radiative trident PDG ID", metavar="isRadPDG", default=622)
base.parser.add_argument("--pSmearingSeed", type=int, dest="pSmearingSeed",
                         help="Set job dependent momentum smearing seed", metavar="pSmearingSeed", default=42)

options = base.parser.parse_args()

# Use the input file to set the output file name
infile = options.inFilename
outfile = options.outFilename

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
vtxana = HpstrConf.Processor('vtxana', 'NewVertexAnaProcessor')
#mcana = HpstrConf.Processor('mcpartana', 'MCAnaProcessor')
###############################
#   Processor Configuration   #
###############################
#RecoHitAna
vtxana.parameters["debug"] = 0
vtxana.parameters["anaName"] = "vtxana"
vtxana.parameters["tsColl"] = "TSBank"
vtxana.parameters["vtxColl"] = "UnconstrainedV0Vertices_KF"
vtxana.parameters["mcColl"] = "MCParticle"
vtxana.parameters["hitColl"] = "SiClustersOnTrackOnPartOnUVtx"
vtxana.parameters["analysis"] = "vertex"
vtxana.parameters["vtxSelectionjson"] = os.environ['HPSTR_BASE']+"/analysis/selections/simps/vertexSelection_2016_simp_preselection.json"
#vtxana.parameters["vtxSelectionjson"] = os.environ['HPSTR_BASE']+"/analysis/selections/simps/vertexSelection_2016_simp_nocuts.json"
vtxana.parameters["histoCfg"] = os.environ['HPSTR_BASE']+"/analysis/plotconfigs/tracking/simps/vtxAnalysis_2016_simp_reach_light.json"
#vtxana.parameters["histoCfg"] = os.environ['HPSTR_BASE']+"/analysis/plotconfigs/tracking/simps/vtxAnalysis_2016_simp_reach.json"
vtxana.parameters["mcHistoCfg"] = os.environ['HPSTR_BASE']+'/analysis/plotconfigs/mc/basicMC.json'
#####
vtxana.parameters["beamE"] = base.beamE[str(options.year)]
vtxana.parameters["isData"] = options.isData
vtxana.parameters["isRadPDG"] = options.isRadPDG
vtxana.parameters["makeFlatTuple"] = options.makeFlatTuple
vtxana.parameters["beamPosCfg"] = ""
if options.isData and options.year == 2016:
    vtxana.parameters["v0ProjectionFitsCfg"] = os.environ['HPSTR_BASE']+'/analysis/data/v0_projection_2016_config.json'
    vtxana.parameters["trackBiasCfg"] = os.environ['HPSTR_BASE']+'/analysis/data/track_bias_corrections_data_2016.json'
elif not options.isData and options.year == 2016:
    print('Running MC')
    vtxana.parameters["trackBiasCfg"] = os.environ['HPSTR_BASE']+'/analysis/data/track_bias_corrections_tritrig_2016.json'
    vtxana.parameters["pSmearingFile"] = os.environ['HPSTR_BASE']+'/utils/data/fee_smearing/smearingFile_2016_all_20240620.root'
    vtxana.parameters["pSmearingSeed"] = options.pSmearingSeed
    vtxana.parameters["v0ProjectionFitsCfg"] = os.environ['HPSTR_BASE']+'/analysis/data/v0_projection_2016_mc_config.json'

CalTimeOffset = -999

if (options.isData == 1):
    CalTimeOffset = 56.
    print("Running on data file: Setting CalTimeOffset %d" % CalTimeOffset)

elif (options.isData == 0):
    CalTimeOffset = 42.4
    print("Running on MC file: Setting CalTimeOffset %d" % CalTimeOffset)
else:
    print("Specify which type of ntuple you are running on: -t 1 [for Data] / -t 0 [for MC]")

vtxana.parameters["CalTimeOffset"] = CalTimeOffset
#Region definitions
RegionPath = os.environ['HPSTR_BASE']+"/analysis/selections/simps/"
if (options.year == 2016):
    RegionDefinitions =  [RegionPath+'Tight_2016_simp_reach_CR.json',
                              RegionPath+'Tight_2016_simp_SR_analysis.json',
                              RegionPath+'Tight_nocuts.json',
                              RegionPath+'Tight_L1L1_nvtx1.json']
    if(options.isData != 1):
        RegionDefinitions.extend([RegionPath+'radMatchTight_2016_simp_reach_CR.json',
            RegionPath+'radMatchTight_2016_simp_SR_analysis.json',
            RegionPath+'radMatchTight_nocuts.json',
            RegionPath+'radMatchTight_L1L1_nvtx1.json']
        )

    vtxana.parameters["regionDefinitions"] = RegionDefinitions

#MCParticleAna
#mcana.parameters["debug"] = 0
#mcana.parameters["anaName"] = "mcAna"
#mcana.parameters["partColl"] = "MCParticle"
#mcana.parameters["trkrHitColl"] = "TrackerHits"
#mcana.parameters["ecalHitColl"] = "EcalHits"
#mcana.parameters["histCfg"] = os.environ['HPSTR_BASE']+'/analysis/plotconfigs/mc/basicMC.json'

# Sequence which the processors will run.
p.sequence = [vtxana]  # ,mcana]

p.skip_events = options.skip_events
p.max_events = options.nevents

p.input_files = infile
p.output_files = [outfile]

p.printProcess()
