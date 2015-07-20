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

#def dataindex_h(radius,angle):
    #return datain[np.logical_and(r == radius, theta == angle)]

def get_coordinate(coord):
    return (resolution/xmax)*coord

matplotlib.rcdefaults()
matplotlib.rc('font',**{'family':'serif','size':16})
os.environ['PATH'] = os.environ['PATH'] + ':/opt/local/bin'
os.environ['PATH']

ymin = 0.#*1477000.
ymax = 100.#*1477000.
xmin = 0.#*1477000.
xmax = 100.#*1477000.
resolution = 100.
savefilename = '/Users/Anton/Desktop/Data/Binaries/hydro_particle_5_50.npy'
onlyplot = False


def plotarrows():
    COORDINATES = np.load(savefilename)
    plt.figure()
    plt.rc('text', usetex=True)
    ypmin = 0
    ypmax = 100
    xpmin = 0
    xpmax = 70
        
    plt.xlim(xmin=xpmin, xmax=xpmax)
    plt.ylim(ymin=ypmin, ymax=ypmax)
    
    for i in range(0,len(COORDINATES)-1):
        if i%2==0:
            plt.arrow(COORDINATES[i,0], COORDINATES[i,1], (COORDINATES[i+1,0] - COORDINATES[i,0]), (COORDINATES[i+1,1] - COORDINATES[i,1]), fc="k", ec="k", head_width=1.5, head_length=1)
        else:
            plt.arrow(COORDINATES[i,0], COORDINATES[i,1], (COORDINATES[i+1,0] - COORDINATES[i,0]), (COORDINATES[i+1,1] - COORDINATES[i,1]))
    
    plt.gca().set_aspect('equal')    
    plt.title('Particle trail')
    plt.xlabel('$r/r_g$')
    plt.ylabel('$r/r_g$')
    plt.show()
    
def iterate():  
    
    #grid_y, grid_x = np.mgrid[ymin:ymax:100j,xmin:xmax:100j]
  
    starting_index    = 1000
    last_index        = 1310 
    TIME_INTERVAL     = 10.
    START_x           = get_coordinate(7)
    START_y           = get_coordinate(50)
    COORDINATES       = np.array([[START_x, START_y]])
    
    coord_x = START_x
    coord_y = START_y
    
    for fileindex in range(starting_index,last_index):
        
        grid_y, grid_x = np.mgrid[coord_y:coord_y:1j,coord_x:coord_x:1j]

        
        if fileindex < 1000:
            filename = '0' + str(fileindex)
        else:
            filename = str(fileindex)
        
        
        filein = open('/Users/Anton/Desktop/Data/hd300a0/hd300a0_rel/sim' + filename + '.dat')
        datain = np.loadtxt(filein,input_dtype.dtype_newtonian())
        filein.close()
        a         = 0
        r         = datain['r']#*1477000.
        theta     = datain['theta']
        u_radial  = datain['u_1']/datain['u_t']#*3.e10
        u_theta   = datain['u_2']/datain['u_t']#*3.e10
        
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
            
            #DELTA_x = grid_ux[get_coordinate(coord_y), get_coordinate(coord_x)]*TIME_INTERVAL
            #DELTA_y = grid_uy[get_coordinate(coord_y), get_coordinate(coord_x)]*TIME_INTERVAL
            
            coord_x = coord_x + DELTA_x
            coord_y = coord_y + DELTA_y
            
            if (np.logical_or(coord_x > xmax, coord_y > ymax)):
                print('The particle have escaped the bounded region. Terminating loop at index ' + str(fileindex))
                break

            COORDINATES = np.append(COORDINATES, np.array([[coord_x, coord_y]]), axis=0)
        except IndexError:
            np.save(savefilename, COORDINATES)
            print('IndexError: The particle have escaped the bounded region. Terminating loop at index ' + str(fileindex))
            break
            
        print('Index:', fileindex, COORDINATES)
        #u_radial_START = grid_u_radial()
        #plt.contourf(grid_x, grid_y, grid_u_radial)
        #
        #s = plt.streamplot(grid_x, grid_y, grid_ux, grid_uy, 
        #        density=1, color='#BBBBBB', 
        #        arrowsize=2.5, arrowstyle='->',
        #        minlength=.1)
        #
        #plt.gca().set_aspect('equal')    
        #plt.show()
    np.save(savefilename, COORDINATES)
    
if onlyplot == True:
    plotarrows()
else:
    iterate()
    plotarrows()








