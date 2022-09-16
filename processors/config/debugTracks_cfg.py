import HpstrConf
import os
import sys
import baseConfig as base

base.parser.add_argument("-w", "--tracking", type=str, dest="tracking",
                          help="Which tracking to use to make plots", metavar="tracking", default="KF")
options = base.parser.parse_args()


# Use the input file to set the output file name
#inFilename  = sys.argv[1].strip()
#outFilename = '%s_anaTrks.root' % inFilename[:-5]
# Use the input file to set the output file name
inFilename = options.inFilename[0]
outFilename = options.outFilename

print('Input file:  %s' % inFilename)
print('Output file: %s' % outFilename)

p = HpstrConf.Process()

p.run_mode = 1
#p.max_events = 1000

# Library containing processors
p.add_library("libprocessors")

###############################
#          Processors         #
###############################
anaTrks = HpstrConf.Processor('anaTrks', 'DebugProcessor')

###############################
#   Processor Configuration   #
###############################
anaTrks.parameters["debug"] = 0 
anaTrks.parameters["isData"] = options.isData
anaTrks.parameters["ecalColl"] = "RecoEcalClusters_GBL"
#anaTrks.parameters["ecalColl"] = "RecoEcalClusters"
anaTrks.parameters["histCfg"] = os.environ['HPSTR_BASE']+'/analysis/plotconfigs/tracking/debug_basicTracking.json'
anaTrks.parameters["selectionjson"] = os.environ['HPSTR_BASE']+'/analysis/selections/trackHit/trackHitAna.json'

CalTimeOffset=-999

if (options.isData==1):
    CalTimeOffset=55.
    print("Running on data file: Setting CalTimeOffset %d"  % CalTimeOffset)

elif (options.isData==0):
    CalTimeOffset=43.
    print("Running on MC file: Setting CalTimeOffset %d"  % CalTimeOffset)
else:
    print("Specify which type of ntuple you are running on: -t 1 [for Data] / -t 0 [for MC]")

anaTrks.parameters["CalTimeOffset"]=CalTimeOffset

if options.tracking == "KF":
    anaTrks.parameters["trkCollName"] = 'KalmanFullTracks'
    anaTrks.parameters["vtxCollName"] = 'UnconstrainedV0Vertices_KF'
elif options.tracking == "GBL":
    anaTrks.parameters["trkCollName"] = 'GBLTracks'
    anaTrks.parameters["vtxCollName"] = 'UnconstrainedV0Vertices'
else:
    print("ERROR! INCORRECT TRACK TYPE SPECIFIED")
anaTrks.parameters["regionDefinitions"] = ['']
# Sequence which the processors will run.
p.sequence = [anaTrks]

p.input_files=[inFilename]
p.output_files = [outFilename]

p.printProcess()
