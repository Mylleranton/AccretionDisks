import matplotlib.pyplot as plt
import matplotlib.lines as mline
from matplotlib.collections import LineCollection
import matplotlib
import os
import numpy as np
from christoffel_symbols import christoffel as ch
import force_vectors_all_forces as fv 


filebase = '/Users/Anton/Desktop/Data/Binaries/hydro_particle'
files = [filebase + '_1000_15_20_dt10.npy',
        filebase + '_1000_10_20_dt10.npy',
        filebase + '_1000_40_30_dt10.npy']

COORD_array = np.empty((1,len(files)))
for i in range(0, len(files)):
    try:
        mLoaded = np.load(files[i])
        COORD_array[i] = mLoaded
    except (IndexError, IOError):
        print('Error. Terminating...')
        exit()
    

plt.figure()
plt.rc('text', usetex=True)
ypmin = 0
ypmax = 100
xpmin = 0
xpmax = 100
    
plt.xlim(xmin=xpmin, xmax=xpmax)
plt.ylim(ymin=ypmin, ymax=ypmax)

for i in range(0,len(COORD_array)):
    mCOORD = COORD_array[i]
    tmp_array = mCOORD[:,0:2].reshape(-1,1,2)
    segments = np.concatenate([tmp_array[:-1], tmp_array[1:]], axis=1)
        
    lw = 0.1+mCOORD[:,2]*200.
    lc = LineCollection(segments, linewidths=lw, colors=['grey'])
    plt.gca().add_collection(lc)
    
    plt.gca().set_aspect('equal')    
    plt.xlabel('$x/r_g$')
    plt.ylabel('$z/r_g$')
