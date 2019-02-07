import sys
sys.path.append("/home/nvaytet/work/code/mantid/builds/py2/bin")
from mantid.simpleapi import *


ws = LoadEmptyInstrument(Filename="V20_geometry.nxs")
inst = ws.getInstrument()
di = ws.detectorInfo()
ci = ws.componentInfo()
print("Modified workspace {0} has instrument: {1}".format(ws.getName(), inst.getName()))
print("Instrument {0} has {1} components, including {2} detectors".format(inst.getName(), ci.size(), di.size()))

ncomps = ci.size()
# detIds = di.detectorIDs()

min_det_x = 1000
max_det_x = 0
min_det_y = 1000
max_det_y = 0
min_det_z = 1000
max_det_z = 0


# Find min and max limits of 2Theta
print("========================")
print("Non-detector components:")
for i in range(ncomps):
    if not ci.isDetector(i):
        print(ci.name(i) + ": ",ci.position(i))
    else:
        min_det_x = min(min_det_x,ci.position(i)[0])
        max_det_x = max(max_det_x,ci.position(i)[0])
        min_det_y = min(min_det_y,ci.position(i)[1])
        max_det_y = max(max_det_y,ci.position(i)[1])
        min_det_z = min(min_det_z,ci.position(i)[2])
        max_det_z = max(max_det_z,ci.position(i)[2])

print("========================")
print("Group of detectors min/max positions:")
print("  x: {} ; {}".format(min_det_x,max_det_x))
print("  y: {} ; {}".format(min_det_y,max_det_y))
print("  z: {} ; {}".format(min_det_z,max_det_z))

print("========================")
print("Source position: {}".format(ci.sourcePosition()))
print("Sample position: {}".format(ci.samplePosition()))
print("========================")
