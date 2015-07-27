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

fig = plt.figure()
plt.rc('text', usetex=True)
ypmin = 0
ypmax = 100
xpmin = 0
xpmax = 100
    
plt.xlim(xmin=xpmin, xmax=xpmax)
plt.ylim(ymin=ypmin, ymax=ypmax)

for i in range(0,len(loadedFiles)):
    mCOORD = loadedFiles[i]
    mCOORD = mCOORD[np.logical_not(mCOORD[:,1]>85.)]
    
    tmp_array = mCOORD[:,0:2].reshape(-1,1,2)
    segments = np.concatenate([tmp_array[:-1], tmp_array[1:]], axis=1)
        
    lw = 0.1+mCOORD[:,2]*200.
    clist = ['red', 'orange', 'gold']
    ls = ['solid', 'dashed', 'dotted']
    norm = matplotlib.colors.Normalize(1.e-4,1.2e-1)
    #lc = LineCollection(segments, linewidths=3, colors=clist[i], linestyles=ls[i])
    lc = LineCollection(segments, linewidths=5, cmap='bwr', norm=norm)
    lc.set_array(mCOORD[:,2])
    plt.gca().add_collection(lc)
    
    plt.arrow(mCOORD[-1,0],mCOORD[-1,1],mCOORD[-1,0]-mCOORD[-5,0],mCOORD[-1,1]-mCOORD[-5,1], 
            width=0.5, head_width=3, overhang=0.3,length_includes_head=False, color='black', zorder=10)
    #string_number = r'\#' + str(i+1)
    #plt.annotate(s=string_number, xy=(mCOORD[0,0]+3*(-1)**(i+1),mCOORD[0,1]-7) )
    
    plt.gca().set_aspect('equal')    
    plt.xlabel('$x/r_g$')
    #plt.ylabel('$z/r_g$')

filein = open('/Volumes/Seagate/4Anton/hd300a0/dt100/simavg0070-0134_rel.dat','rb')
datain = np.loadtxt(filein,dt.dtype_rel_hydro())
filein.close()
grid_y, grid_x = np.mgrid[ypmin:ypmax:200j,xpmin:xpmax:200j]
points = (datain['r']*np.sin(datain['theta']),datain['r']*np.cos(datain['theta']))
grid_rho = griddata(points, datain['rho'], (grid_x, grid_y), method='linear')


ax = fig.gca()

cax, args = matplotlib.colorbar.make_axes(ax, anchor=(-1.4, 0.95), shrink=0.55)
colorbar = plt.colorbar(lc, ticks=np.linspace(0,0.5,11), cax=cax)
plt.sca(ax)


c = plt.contourf(grid_x, grid_y, np.log10(grid_rho*6.17e15), extend='both', levels=np.linspace(-8,-4,61), cmap='gray', alpha=0.99)
colorbar2 = plt.colorbar(c, ticks=np.linspace(-8,-4,5).astype(int))


p1 = mline.Line2D([], [], color='blue', label='Gravity', linewidth=2.5)
p2 = mline.Line2D([], [], color='green', label='Thermal', linewidth=2.5)
p3 = mline.Line2D([], [], color='darkmagenta', label='Magnetic', linewidth=2.5)
p4 = mline.Line2D([], [], color='red', label='Centrifugal', linewidth=2.5)
p5 = mline.Line2D([], [], color='cyan', label='Rel. correction', linewidth=2.5)
p7 = mline.Line2D([], [], color='black', label='Total', linewidth=2.5)

legend = plt.legend(handles=[p1,p2,p3,p4,p5,p7], fontsize='small',loc='lower right', bbox_to_anchor=(1, 0.01))
legend.draw_frame(False)

plt.show()
plt.savefig('/Users/Anton/Dropbox/Aleksander/Figures/simavg0070-0134/particles/three_particles', bbox_inches='tight') 
