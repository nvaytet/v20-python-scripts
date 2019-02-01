import sys
sys.path.append("/home/nvaytet/work/code/mantid/builds/py2/bin")
from mantid.simpleapi import *

# theFile = "test.nxs"
# theFile = "V20_ESSIntegration_2018-12-11_1743_agg.nxs"
# theFile = "/home/nvaytet/work/code/testing-framework/python/V20_ESSIntegration_2018-12-13_0942_agg.nxs"
theFile = "V20_ESSIntegration_2018-12-12_1209_agg.nxs"



ws1 = CreateSampleWorkspace();
mon1 = LoadInstrument(ws1, FileName=theFile, RewriteSpectraMap=True)
inst1 = ws1.getInstrument()
di1 = ws1.detectorInfo()
ci1 = ws1.componentInfo()
print("Modified workspace {0} has instrument: {1}".format(ws1.getName(), inst1.getName()))
print("Instrument {0} has {1} components, including {2} monitors and {3} detectors".format(inst1.getName(), ci1.size(), len(mon1), di1.size()))

# ws2 = LoadEventNexus(FileName=theFile, LoadLogs=True)
