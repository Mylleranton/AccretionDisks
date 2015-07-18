import matplotlib.pyplot as plt
import matplotlib
import os
import numpy as np
from scipy.interpolate import griddata 

#Read, plot density (interapolate), vector field, streamlines
# Magnetic pressure (bsq) and direction of B-field as streamlines

#Plot preamble
matplotlib.rcdefaults()
matplotlib.rc('font',**{'family':'serif','size':16})
os.environ['PATH'] = os.environ['PATH'] + ':/opt/local/bin'
os.environ['PATH']

#File operations
filein = open('/Users/Anton/Desktop/Data/simavg0070-0134.dat','rb')
fileout = '/Users/Anton/Desktop/Data/PP_001.png'

input_dtype = np.dtype([
    ("ix", int),
    ("iy", int),
    ("iz", int),
    ("r", float),
    ("theta", float),
    ("phi", float),
    ("rho", float),
    ("temp", float),
    ("u_t", float),
    ("u_1", float),
    ("u_2", float),
    ("u_3", float),
    ("volume", float),
    ("bsq", float),
    ("b_1", float),
    ("b_2", float),
    ("b_3", float)
    ])

datain = np.loadtxt(filein,input_dtype)
filein.close()

#Interpolate data
ymin = -50
ymax = 50
xmin = 0
xmax = 150

grid_y, grid_x = np.mgrid[ymin:ymax:150j,xmin:xmax:150j]

#Convert coordinates from spherical to cartesian
points = (datain['r']*np.sin(datain['theta']),datain['r']*np.cos(datain['theta']))

ux = (-datain['u_2']*np.cos(datain['theta'])*datain['r'] + datain['u_1']*np.sin(datain['theta']))
uy = (-datain['u_2']*np.sin(datain['theta'])*datain['r'] + datain['u_1']*np.cos(datain['theta']))

grid_rho = griddata(points, datain['rho'], (grid_x, grid_y), method='linear')
grid_ut  = griddata(points, datain['u_t'], (grid_x, grid_y), method='linear')
grid_ux  = griddata(points, ux, (grid_x, grid_y), method='linear')
grid_uy  = griddata(points, uy, (grid_x, grid_y), method='linear')

#proper_v_1 = np.divide(grid_ux, grid_ut)
#proper_v_2 = np.divide(grid_uy, grid_ut)

# Plot the interpolated data
plt.close('all')
plt.figure()
plt.rc('text', usetex=True)

c = plt.contourf(grid_x, grid_y, grid_rho, 150, extend='both' 
    #levels=np.linspace(-16,-8,60)
    )
plt.colorbar(c)
plt.title('Density distribution and particle velocity [CGS]')
plt.xlabel('$r / r_g$')
plt.ylabel('$r / r_g$')

s = plt.streamplot(grid_x, grid_y, grid_ux, grid_uy, 
        density=1, color='#FFFFFF', 
        arrowsize=2.5, arrowstyle='->',
        minlength=.1)

plt.savefig(fileout)
plt.show()