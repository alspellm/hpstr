import HpstrConf
import sys

import baseConfig as base
from baseConfig import bfield

options = base.parser.parse_args()

# Use the input file to set the output file name
lcio_file = options.inFilename
root_file = options.outFilename

print('LCIO file: %s' % lcio_file)
print('Root file: %s' % root_file)

p = HpstrConf.Process()

#p.max_events = 1000
p.run_mode = 0

# Library containing processors
p.add_library("libprocessors")

###############################
#          Processors         #
###############################
header  = HpstrConf.Processor('header', 'EventProcessor')
track   = HpstrConf.Processor('track', 'TrackingProcessor')
trackgbl = HpstrConf.Processor('trackgbl', 'TrackingProcessor')
trackrefitgbl = HpstrConf.Processor('trackrefitgbl', 'TrackingProcessor')
svthits = HpstrConf.Processor('svthits', 'Tracker3DHitProcessor')
rawsvt  = HpstrConf.Processor('rawsvt', 'SvtRawDataProcessor')
ecal    = HpstrConf.Processor('ecal', 'ECalDataProcessor')
vtx     = HpstrConf.Processor('vtx', 'VertexProcessor')
vtxgbl   = HpstrConf.Processor('vtxgbl', 'VertexProcessor')
cvtxgbl   = HpstrConf.Processor('cvtxgbl', 'VertexProcessor')
mcpart  = HpstrConf.Processor('mcpart', 'MCParticleProcessor')

###############################
#   Processor Configuration   #
###############################
#Event
header.parameters["debug"] = 0
header.parameters["headCollRoot"] = "EventHeader"
header.parameters["trigCollLcio"] = "TriggerBank"
header.parameters["rfCollLcio"]   = "RFHits"
header.parameters["vtpCollLcio"]  = "VTPBank"
header.parameters["vtpCollRoot"]  = "VTPBank"
header.parameters["tsCollLcio"]   = "TSBank"
header.parameters["tsCollRoot"]   = "TSBank"

#SvtRawData
rawsvt.parameters["debug"] = 0
rawsvt.parameters["hitCollLcio"]    = 'SVTRawTrackerHits'
rawsvt.parameters["hitfitCollLcio"] = 'SVTFittedRawTrackerHits'
rawsvt.parameters["hitCollRoot"]    = 'SVTRawTrackerHits'

#Tracker3DHits
svthits.parameters["debug"] = 0
svthits.parameters["hitCollLcio"]    = 'RotatedHelicalTrackHits'
svthits.parameters["hitCollRoot"]    = 'RotatedHelicalTrackHits'


#Tracking
track.parameters["debug"] = 0
track.parameters["trkCollLcio"] = 'KalmanFullTracks'
track.parameters["trkCollRoot"] = 'KalmanFullTracks'
track.parameters["kinkRelCollLcio"] = ''
track.parameters["trkRelCollLcio"] = 'KFTrackDataRelations'
track.parameters["trkhitCollRoot"] = 'SiClustersOnTrack'
track.parameters["hitFitsCollLcio"] = 'SVTFittedRawTrackerHits'
track.parameters["rawhitCollRoot"] = ''

#Only for detail studies
#LT uncomment
track.parameters["rawhitCollRoot"] = ''#'SCTRawHitsOnTrack_KF'

#LT uncommented
if (not options.isData):
    track.parameters["truthTrackCollLcio"] = 'KalmanFullTracksToTruthTrackRelations'
    track.parameters["truthTrackCollRoot"] = 'KalmanFullTracksTruth'
    track.parameters["trackTruthInfoCollLcio"] = 'KalmanFullTracksToTruthInfoRelations'
RegionPath = os.environ['HPSTR_BASE']+"/analysis/selections/trackHit/"
anaTrks.parameters["regionDefinitions"] = [RegionPath+'hc15_1111.json',
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
                                           RegionPath+'noHCreq.json'
                                          ] 

#LT check if we need the b field or not -- version of HPS java
# for Jess's files need to give it b-field

track.parameters["bfield"] = bfield[str(options.year)]

trackgbl.parameters["debug"] = 0
trackgbl.parameters["trkCollLcio"] = 'GBLTracks'
trackgbl.parameters["trkCollRoot"] = 'GBLTracks'
trackgbl.parameters["kinkRelCollLcio"] = 'GBLKinkDataRelations'
trackgbl.parameters["trkRelCollLcio"] = 'TrackDataRelations'
trackgbl.parameters["trkhitCollRoot"] = 'RotatedHelicalOnTrackHits'
trackgbl.parameters["hitFitsCollLcio"] = 'SVTFittedRawTrackerHits'

#Only for detail studies
#LT uncomment
trackgbl.parameters["rawhitCollRoot"] = ''#'SVTRawHitsOnTrack'

#LT uncommented
if (not options.isData):
    trackgbl.parameters["truthTrackCollLcio"] = 'GBLTracksToTruthTrackRelations'
    trackgbl.parameters["truthTrackCollRoot"] = 'GBLTracksTruth'

#LT check if we need the b field or not -- version of HPS java
trackgbl.parameters["bfield"] = bfield[str(options.year)]

#ECalData
ecal.parameters["debug"] = 0
ecal.parameters["hitCollLcio"] = 'EcalCalHits'
ecal.parameters["hitCollRoot"] = ''#'RecoEcalHits'
ecal.parameters["clusCollLcio"] = "EcalClustersCorr"
ecal.parameters["clusCollRoot"] = "RecoEcalClusters"

#Vertex
vtx.parameters["debug"] = 0
vtx.parameters["vtxCollLcio"]    = 'UnconstrainedV0Vertices_KF'
vtx.parameters["vtxCollRoot"]    = 'UnconstrainedV0Vertices_KF'
vtx.parameters["partCollRoot"]   = 'ParticlesOnVertices_KF'
vtx.parameters["kinkRelCollLcio"] = ''
vtx.parameters["trkRelCollLcio"] = 'KFTrackDataRelations'


vtxgbl.parameters["debug"] = 0
vtxgbl.parameters["vtxCollLcio"]     = 'UnconstrainedV0Vertices'
vtxgbl.parameters["vtxCollRoot"]     = 'UnconstrainedV0Vertices'
vtxgbl.parameters["partCollRoot"]    = 'ParticlesOnVertices'
vtxgbl.parameters["kinkRelCollLcio"] = 'GBLKinkDataRelations'
vtxgbl.parameters["trkRelCollLcio"]  = 'TrackDataRelations'


cvtxgbl.parameters["debug"] = 0
cvtxgbl.parameters["vtxCollLcio"]     = 'TargetConstrainedV0Vertices'
cvtxgbl.parameters["vtxCollRoot"]     = 'TargetConstrainedV0Vertices'
cvtxgbl.parameters["partCollRoot"]    = 'ParticlesOnVertices'
cvtxgbl.parameters["kinkRelCollLcio"] = 'GBLKinkDataRelations'
cvtxgbl.parameters["trkRelCollLcio"]  = 'TrackDataRelations'


#MCParticle
mcpart.parameters["debug"] = 0
mcpart.parameters["mcPartCollLcio"] = 'MCParticle'
mcpart.parameters["mcPartCollRoot"] = 'MCParticle'

# Sequence which the processors will run.
if (not options.isData):
    p.sequence = [header,track]
else:
    p.sequence = [header,track]

if (options.nevents > -1 ):
    p.max_events = options.nevents


p.input_files = lcio_file
p.output_files = [root_file]



p.printProcess()
