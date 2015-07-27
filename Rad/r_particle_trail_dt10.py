import matplotlib.pyplot as plt
import matplotlib.lines as mline
import matplotlib
import os
import numpy as np
from scipy.interpolate import griddata 
from christoffel_symbols import christoffel as ch
import STATIC_DTYPE as input_dtype

def dataindex(ix, iy):
    return ix + 252*iy

def get_coordinate(coord):
    return (float(resolution)/float(xmax))*float(coord)

matplotlib.rcdefaults()
matplotlib.rc('font',**{'family':'serif','size':16})
os.environ['PATH'] = os.environ['PATH'] + ':/opt/local/bin'
os.environ['PATH']

ymin = 0.
ymax = 100.
xmin = 0.
xmax = 100.
resolution = 100.
gSTART_x = 69
gSTART_y = 70
gSTART_INDEX = 1000
gLAST_INDEX = 1388
loop_index = 0

onlyplot = True
REVERSE = False

savefilename = str(gSTART_INDEX) + '_' + str(gSTART_x)+'_'+str(gSTART_y)+'_dt10'
savebase = '/Users/Anton/Desktop/Data/Binaries/radiative_particle_'
figbase = '/Users/Anton/Dropbox/Aleksander/Figures/simavg0100-0189/particles/particle_'

if REVERSE == True:
    savebase = savebase + 'R_'
    figbase = figbase + 'R_'

def plotarrows():
    COORDINATES = np.load(savebase + savefilename + '.npy')
    plt.figure()
    plt.rc('text', usetex=True)
    ypmin = 0
    ypmax = 100
    xpmin = 0
    xpmax = 70
        
    plt.xlim(xmin=xpmin, xmax=xpmax)
    plt.ylim(ymin=ypmin, ymax=ypmax)
    
    for i in range(0,len(COORDINATES)-1):
        lw = 0.1+COORDINATES[i,2]*100.
        if i%21==0:
            plt.arrow(COORDINATES[i,0], COORDINATES[i,1], (COORDINATES[i+1,0] - COORDINATES[i,0]), (COORDINATES[i+1,1] - COORDINATES[i,1]), fc="k", ec="k", head_width=1.5, head_length=1, linewidth=lw)
        else:
            plt.arrow(COORDINATES[i,0], COORDINATES[i,1], (COORDINATES[i+1,0] - COORDINATES[i,0]), (COORDINATES[i+1,1] - COORDINATES[i,1]), linewidth=lw)
    
    plt.gca().set_aspect('equal')    
    plt.title('Trajectory over time: ' + str(len(COORDINATES)*10))
    plt.xlabel('$x/r_g$')
    plt.ylabel('$z/r_g$')

    plt.show()
    plt.savefig(figbase+savefilename, bbox_inches='tight') 

    
def iterate():    
    starting_index    = gSTART_INDEX
    last_index        = gLAST_INDEX 
    TIME_INTERVAL     = 10.
    START_x           = get_coordinate(gSTART_x)
    START_y           = get_coordinate(gSTART_y)
    COORDINATES       = np.array([[START_x, START_y, 0]])
    
    coord_x = START_x
    coord_y = START_y
    
    for fileindex in range(starting_index,last_index):
        grid_y, grid_x = np.mgrid[coord_y:coord_y:1j,coord_x:coord_x:1j]

        
        if fileindex < 1000:
            filename = '0' + str(fileindex)
        else:
            filename = str(fileindex)
        
        
        #filein = open('/Users/Anton/Desktop/Data/hd300a0/hd300a0_rel/sim' + filename + '.dat')
        filein = open('/Volumes/Seagate/4Anton/d300a0/dt10/simext' + filename + '.dat')
        datain = np.loadtxt(filein,input_dtype.dtype_rel_radiation())
        filein.close()
        a         = 0
        r         = datain['r']
        theta     = datain['theta']
        u_radial  = datain['u_1']/datain['u_t']
        u_theta   = datain['u_2']/datain['u_t']
        
        g11 = np.sqrt((np.divide( (np.power(r, 2)+a**(2.)*np.power(np.cos(theta),2)), (np.power(r, 2)-2.*r+a**2) ) ))
        g22 = np.sqrt((np.power(r,2)+a**(2.)*np.power(np.cos(theta),2)))
        
        points = (r*np.sin(theta),r*np.cos(theta))
        ux = (u_theta*np.cos(theta)*g22 + u_radial*np.sin(theta))
        uy = (-u_theta*np.sin(theta)*g22 + u_radial*np.cos(theta))
        
        grid_ux = griddata(points, ux, (grid_x, grid_y), method='linear', fill_value=1e-30)
        grid_uy = griddata(points, uy, (grid_x, grid_y), method='linear', fill_value=1e-30)
        try: 
            DELTA_x = grid_ux[0,0]*TIME_INTERVAL
            DELTA_y = grid_uy[0,0]*TIME_INTERVAL
            coord_x = coord_x + DELTA_x
            coord_y = coord_y + DELTA_y
            
            if (np.logical_or(coord_x > xmax, coord_y > ymax)):
                print('The particle have escaped the bounded region. Terminating loop at index ' + str(fileindex))
                break
            elif (np.logical_and(coord_x < 2, coord_y < 2)):
                print('The particle have been engulfed by the black hole. Terminating loop at index ' + str(fileindex))
                break
                
            COORDINATES = np.append(COORDINATES, np.array([[coord_x, coord_y, np.sqrt(grid_ux[0,0]**(2.) + grid_uy[0,0]**(2.))]]), axis=0)
        except IndexError:
            np.save(savefilename, COORDINATES)
            print('IndexError: The particle have escaped the bounded region. Terminating loop at index ' + str(fileindex))
            break
            
        print('Index:', fileindex, COORDINATES)
    np.save(savebase + savefilename + '.npy', COORDINATES)
    
def iterate_reverse():
    starting_index    = gSTART_INDEX
    last_index        = gLAST_INDEX 
    TIME_INTERVAL     = 10.
    START_x           = get_coordinate(gSTART_x)
    START_y           = get_coordinate(gSTART_y)
    COORDINATES       = np.array([[START_x, START_y, 0]])
    
    coord_x = START_x
    coord_y = START_y
    
    for loop_index in range(0,300):         
        fileindex = starting_index - loop_index
        
        grid_y, grid_x = np.mgrid[coord_y:coord_y:1j,coord_x:coord_x:1j]

        
        if fileindex < 1000:
            filename = '0' + str(fileindex)
        else:
            filename = str(fileindex)
        
        
        #filein = open('/Users/Anton/Desktop/Data/hd300a0/hd300a0_rel/sim' + filename + '.dat')
        filein = open('/Volumes/Seagate/4Anton/d300a0/dt10/simext' + filename + '.dat')
        datain = np.loadtxt(filein,input_dtype.dtype_rel_radiation())
        filein.close()
        a         = 0
        r         = datain['r']
        theta     = datain['theta']
        u_radial  = datain['u_1']/datain['u_t']
        u_theta   = datain['u_2']/datain['u_t']
        
        g11 = np.sqrt((np.divide( (np.power(r, 2)+a**(2.)*np.power(np.cos(theta),2)), (np.power(r, 2)-2.*r+a**2) ) ))
        g22 = np.sqrt((np.power(r,2)+a**(2.)*np.power(np.cos(theta),2)))
        
        points = (r*np.sin(theta),r*np.cos(theta))
        ux = (u_theta*np.cos(theta)*g22 + u_radial*np.sin(theta))
        uy = (-u_theta*np.sin(theta)*g22 + u_radial*np.cos(theta))
        
        grid_ux = griddata(points, ux, (grid_x, grid_y), method='linear', fill_value=1e-30)
        grid_uy = griddata(points, uy, (grid_x, grid_y), method='linear', fill_value=1e-30)
        try: 
            DELTA_x = grid_ux[0,0]*TIME_INTERVAL
            DELTA_y = grid_uy[0,0]*TIME_INTERVAL
            coord_x = coord_x - DELTA_x
            coord_y = coord_y - DELTA_y
            
            if (np.logical_or(coord_x > xmax, coord_y > ymax)):
                print('The particle have escaped the bounded region. Terminating loop at index ' + str(fileindex))
                break
            elif (np.logical_and(coord_x < 2, coord_y < 2)):
                print('The particle have been engulfed by the black hole. Terminating loop at index ' + str(fileindex))
                break
                
            COORDINATES = np.append(COORDINATES, np.array([[coord_x, coord_y, np.sqrt(grid_ux[0,0]**(2.) + grid_uy[0,0]**(2.))]]), axis=0)
        except IndexError:
            np.save(savefilename, COORDINATES)
            print('IndexError: The particle have escaped the bounded region. Terminating loop at index ' + str(fileindex))
            break
            
        print('Index:', fileindex, COORDINATES)
    np.save(savebase + savefilename + '.npy', COORDINATES)
    
if onlyplot == True:
    plotarrows()
else:
    if REVERSE == True:    
        iterate_reverse()
        plotarrows()
    else:
        iterate()
        plotarrows()








