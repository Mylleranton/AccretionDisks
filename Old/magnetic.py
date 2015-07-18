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
fileout = '/Users/Anton/Desktop/Data/image_003.png'

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
xmax = 100

grid_y, grid_x = np.mgrid[ymin:ymax:150j,xmin:xmax:150j]

#Convert coordinates from spherical to cartesian
points = (datain['r']*np.sin(datain['theta']),datain['r']*np.cos(datain['theta']))

bx = (-datain['b_2']*np.cos(datain['theta'])*datain['r'] + datain['b_1']*np.sin(datain['theta']))
by = (-datain['b_2']*np.sin(datain['theta'])*datain['r'] + datain['b_1']*np.cos(datain['theta']))

#thermal_pressure = (np.divide(1.38,1.67)*np.multiply(datain['rho'],datain['temp'])*np.power(10.,8)) # Boltzmann*rho*temp/(mu*protonmass)
thermal_pressure = datain['rho']*datain['temp']*1.38e8/1.67 # Boltzmann*rho*temp/(mu*protonmass)
pressure_ratio = np.divide((0.5*datain['bsq']), thermal_pressure)

grid_ratio= griddata(points, pressure_ratio, (grid_x, grid_y), method='linear')
grid_bx  = griddata(points, bx, (grid_x, grid_y), method='linear')
grid_by  = griddata(points, by, (grid_x, grid_y), method='linear')


# Plot the interpolated data
#plt.close('all')
plt.figure()
plt.rc('text', usetex=True)

c = plt.contourf(grid_x, grid_y, np.log10(grid_ratio), levels=np.linspace(-3,2,10), extend='both')
plt.colorbar(c)
plt.title('Magnetic versus thermal pressure and B-field')
plt.xlabel('r / $r_s$')
plt.ylabel('r / $r_s$')

s = plt.streamplot(grid_x, grid_y, grid_bx, grid_by, 
        density=2, color='#dbe8e8', 
        arrowsize=2, arrowstyle='->',
        minlength=0.5)

plt.savefig(fileout)
plt.show()