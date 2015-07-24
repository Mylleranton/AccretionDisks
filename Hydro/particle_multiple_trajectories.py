import matplotlib.pyplot as plt
import matplotlib.lines as mline
from matplotlib.collections import LineCollection
import matplotlib
import os
import numpy as np
from christoffel_symbols import christoffel as ch
import STATIC_DTYPE as dt
from scipy.interpolate import griddata 




filebase = '/Users/Anton/Desktop/Data/Binaries/hydro_particle'
files = [filebase + '_1000_15_20_dt10.npy',
        filebase + '_1000_10_20_dt10.npy',
        filebase + '_1000_69_80_dt10.npy',
        filebase + '_1000_50_90_dt10.npy']

loadedFiles = [np.load(files[0]), np.load(files[1]),np.load(files[2]),np.load(files[3])]

plt.figure()
plt.rc('text', usetex=True)
ypmin = 0
ypmax = 100
xpmin = 0
xpmax = 100
    
plt.xlim(xmin=xpmin, xmax=xpmax)
plt.ylim(ymin=ypmin, ymax=ypmax)

for i in range(0,len(loadedFiles)):
    mCOORD = loadedFiles[i]
    tmp_array = mCOORD[:,0:2].reshape(-1,1,2)
    segments = np.concatenate([tmp_array[:-1], tmp_array[1:]], axis=1)
        
    lw = 0.1+mCOORD[:,2]*200.
    clist = ['red', 'orange', 'gold']
    ls = ['solid', 'dashed', 'dotted']
    norm = matplotlib.colors.Normalize(1.e-4,1.5e-1)
    #lc = LineCollection(segments, linewidths=3, colors=clist[i], linestyles=ls[i])
    lc = LineCollection(segments, linewidths=5, cmap='hot_r', norm=norm)
    lc.set_array(mCOORD[:,2])
    plt.gca().add_collection(lc)
    
    plt.gca().set_aspect('equal')    
    plt.xlabel('$x/r_g$')
    #plt.ylabel('$z/r_g$')

filein = open('/Volumes/Seagate/4Anton/hd300a0/dt100/simavg0070-0134_rel.dat','rb')
datain = np.loadtxt(filein,dt.dtype_rel_hydro())
filein.close()
grid_y, grid_x = np.mgrid[ypmin:ypmax:200j,xpmin:xpmax:200j]
points = (datain['r']*np.sin(datain['theta']),datain['r']*np.cos(datain['theta']))
grid_rho = griddata(points, datain['rho'], (grid_x, grid_y), method='linear')
#colorbar = plt.colorbar(lc)

c = plt.contourf(grid_x, grid_y, np.log10(grid_rho), extend='both', levels=np.linspace(-24,-20,21), cmap='gray', alpha=0.99)


plt.show()
plt.savefig('/Users/Anton/Dropbox/Aleksander/Figures/simavg0070-0134/particles/three_particles', bbox_inches='tight') 
