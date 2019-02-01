import numpy as np
import matplotlib.pyplot as plt

# Load data file
theFile = "V20_ESSIntegration_2018-12-12_1209_agg.txt"

# Get dimensions and extents
with open(theFile, 'r') as f:
    content = f.readline()
    dims = content.split(":")[1].replace("(","").replace(")","").split(',')
    content = f.readline()
    extents = content.split(":")[1].split()
# Convert strings to numbers
nx,ny,nz = np.asarray(dims, dtype=int)
xmin,xmax,ymin,ymax,zmin,zmax = np.asarray(extents, dtype=float)

# Load the whole array
data = np.loadtxt(theFile)
data = data.reshape((nx,ny,nz))

# Make a 2D array (x vs tof) by summing over y
image = np.zeros([nx,nz])
for i in range(nx):
    for k in range(nz):
        image[i,k] = np.sum(data[i,:,k])

# Make figure
fig = plt.figure()
ratio = 1.0
sizex = 10.0
fig.set_size_inches(sizex,ratio*sizex)
# Plot image on top panel
ax1 = fig.add_axes([0, 0.53, 1.0, 0.47])
contf = ax1.imshow(image, origin='lower',extent=[zmin,zmax,xmin,xmax], aspect='auto')
ax1.set_xlabel("Arrival time at detector (microseconds)")
ax1.set_ylabel("Detector x position (cm)")
ax3 = fig.add_axes([1.01, 0.53, 0.03, 0.47])
cb = plt.colorbar(contf,ax=ax1,cax=ax3)
# Plot data summer over all detectors on bottom panel
xsum = np.zeros([nz])
for i in range(nz):
    xsum[i] = np.sum(image[:,i])
dz = (zmax-zmin) / nz
ax2 = fig.add_axes([0, 0.0, 1.0, 0.47])
ax2.plot(np.linspace(zmin+0.5*dz,zmax-0.5*dz,nz),xsum)
ax2.set_xlim([zmin,zmax])
ax2.set_xlabel("Arrival time at detector (microseconds)")
ax2.set_ylabel("Counts (integrated over x)")
# Save figure to file
fig.savefig('image.png',bbox_inches='tight')
