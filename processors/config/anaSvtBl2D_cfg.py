import HpstrConf
import os
import baseConfig as base


base.parser.add_argument('--runNumber', '-r',type=int, dest="runNumber", help="Provide run number.", default="0")
options = base.parser.parse_args()

# Use the input file to set the output file name
lcio_file = options.inFilename[0]
root_file = options.outFilename

print('LCIO file: %s' % lcio_file)
print('Root file: %s' % root_file)

p = HpstrConf.Process()
p.run_mode = 1

#Set files to process
p.input_files=[lcio_file]
p.output_files = [root_file]

# Library containing processors
p.libraries.append("libprocessors.so")

###############################
#          Processors         #
###############################

svtblana = HpstrConf.Processor('svtblana', 'SvtBl2DAnaProcessor')

###############################
#   Processor Configuration   #
###############################
#SvtBl2D
svtblana.parameters["debug"] = 0
svtblana.parameters["removeBaseline"] = 0
svtblana.parameters["baselineFits"] = "filename"
svtblana.parameters["runNumber"] = options.runNumber
svtblana.parameters["rawSvtHitsColl"] = "SVTRawTrackerHits"
svtblana.parameters["histCfg"] = os.environ['HPSTR_BASE']+'/analysis/plotconfigs/svt/Svt2DBl.json'

#p.max_events = 1000

# Sequence which the processors will run.
p.sequence = [svtblana]

p.printProcess()
