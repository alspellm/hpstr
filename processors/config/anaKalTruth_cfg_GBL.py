import HpstrConf
import sys
import os
import baseConfig as base

options = base.parser.parse_args()


# Use the input file to set the output file name
infile = options.inFilename
outfile = options.outFilename

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
trkana    = HpstrConf.Processor('trkana','TruthTrackingAnaProcessor')
trkgblana = HpstrConf.Processor('trkgblana','TruthTrackingAnaProcessor')

trkana.parameters["debug"]  = 0;
trkana.parameters["trkCollName"] = "KalmanFullTracks"
trkana.parameters["histCfg"] = os.environ['HPSTR_BASE'] + '/analysis/plotconfigs/tracking/basicTracking.json'
trkana.parameters["doTruth"] = 1
trkana.parameters["truthHistCfg"] = os.environ['HPSTR_BASE'] + '/analysis/plotconfigs/tracking/truthTrackComparison.json'


trkgblana.parameters["debug"]        = 0;
trkgblana.parameters["trkCollName"]  = "GBLTracks"
trkgblana.parameters["histCfg"]      = os.environ['HPSTR_BASE'] + '/analysis/plotconfigs/tracking/basicTracking.json'
trkgblana.parameters["doTruth"]      = 1
trkgblana.parameters["truthHistCfg"] = os.environ['HPSTR_BASE'] + '/analysis/plotconfigs/tracking/truthTrackComparison.json'
trkgblana.parameters["puritycut"] = 1.0
RegionPath = os.environ['HPSTR_BASE']+"/analysis/selections/trackHit/"
trkgblana.parameters["regionDefinitions"] = [RegionPath+'hc15_1111.json',
                                           RegionPath+'hc14_1110.json',
                                           RegionPath+'hc13_1101.json',
                                           RegionPath+'hc12_1100.json',
                                           RegionPath+'hc11_1011.json',
                                           RegionPath+'hc10_1010.json',
                                           RegionPath+'hc9_1001.json',
                                           RegionPath+'hc8_1000.json',
                                           RegionPath+'hc7_0111.json',
                                           RegionPath+'hc6_0110.json',
                                           RegionPath+'hc5_0101.json',
                                           RegionPath+'hc4_0100.json',
                                           RegionPath+'hc3_0011.json',
                                           RegionPath+'hc2_0010.json',
                                           RegionPath+'hc1_0001.json',
                                           RegionPath+'hc0_0000.json',
                                          ]

p.sequence = [trkgblana]

p.input_files = infile
p.output_files = [outfile]

p.printProcess()
