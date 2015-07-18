import matplotlib.pyplot as plt
import matplotlib.lines as mline
import matplotlib
import os
import numpy as np
from scipy.interpolate import griddata 

#Plot preamble
matplotlib.rcdefaults()
matplotlib.rc('font',**{'family':'serif','size':16})
os.environ['PATH'] = os.environ['PATH'] + ':/opt/local/bin'
os.environ['PATH']

#File operations
index = 700
filename = 'sim0' + str(index)
filein = open('/Users/Anton/Desktop/Data/d300a0/' + filename + '.dat','rb')
fileout = '/Users/Anton/Desktop/Data/Image/' + filename + '.png'

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

# Create meshgrid (in units of gravitational radii, then converted to physical)
ymin = 0
ymax = 100
xmin = 0
xmax = 100

grid_y, grid_x = np.mgrid[ymin:ymax:150j,xmin:xmax:150j]
grid_y = 1477000*grid_y
grid_x = 1477000*grid_x

#Convert coordinates from spherical to cartesian (physical units)
points = (1477000*datain['r']*np.sin(datain['theta']),1477000*datain['r']*np.cos(datain['theta']))

####### Forces ########
## Gravitational force ##
#
# G = 6.674e-8 cm3/g/s2
# M-sun = 1.9891e33 g
F_gravity = 6.674*np.power(10.00, -8)*10*1.9891*np.power(10.00, 33)/(np.square(1477000*datain['r']))

## Thermal force (gas pressure) ##
thermal_pressure = (np.divide(1.38,1.67)*np.multiply(datain['rho'],datain['temp'])*np.power(10.,8)) # Boltzmann*rho*temp/(mu*protonmass)
F_gas = -np.divide(1, datain['rho'])*np.gradient(thermal_pressure)

## Magnetic force (magnetic pressure) ##
#F_mag_vx = (-datain['b_2']*np.cos(datain['theta'])*datain['r'] + datain['b_1']*np.sin(datain['theta']))
#F_mag_vy = (-datain['b_2']*np.sin(datain['theta'])*datain['r'] + datain['b_1']*np.cos(datain['theta']))
F_magnetic = -np.divide(1, datain['rho'])*np.divide(datain['bsq'], 2)

## Centrifugal force ##
v_phi = np.divide(datain['u_3'], datain['u_t'])*np.sqrt(np.square(1477000*datain['r'])*np.square(np.sin(datain['theta'])))
effective_radius = 1477000*datain['r']*np.sin(datain['theta'])
F_centrifugal = np.divide(np.square(v_phi), effective_radius)


##### Interpolate forces and vector fields #####
grid_F_gravity     = griddata(points, F_gravity, (grid_x, grid_y), method='linear')
grid_F_gas         = griddata(points, F_gas, (grid_x, grid_y), method='linear')
grid_F_magnetic    = griddata(points, F_magnetic, (grid_x, grid_y), method='linear')
grid_F_centrifugal = griddata(points, F_centrifugal, (grid_x, grid_y), method='linear')
grid_rho           = griddata(points, datain['rho'], (grid_x, grid_y), method='linear')

########### ------------------------------------------------------------------------------------ ###########
########### ------------------------------------------------------------------------------------ ###########
########### ------------------------------------------------------------------------------------ ###########


# Wind = (40,70), Jet = (15,80), Disk = (80,20)
x         = np.array([10,50,40])
X_jet     = np.array([x[0],x[0],x[0],x[0],x[0]])
X_wind    = np.array([x[1],x[1],x[1],x[1],x[1]])
X_disk    = np.array([x[2],x[2],x[2],x[2],x[2]])
y         = np.array([40,75,15])
Y_jet     = np.array([y[0],y[0],y[0],y[0],y[0]])
Y_wind    = np.array([y[1],y[1],y[1],y[1],y[1]])
Y_disk    = np.array([y[2],y[2],y[2],y[2],y[2]])

#### Returns the logarithm of the value while preserving the sign.
def scale(value):
    tmp_array = value
    #for i in range(0, len(value)):
    #    if value[i] >= 1:       #+1
    #        tmp_array[i] = np.log10(value[i])
    #        continue
    #    elif value[i] > 0:
    #        tmp_array[i] = -np.log10(np.absolute(value[i]))
    #        continue
    #    elif value[i] < 0:    #-
    #        tmp_array[i] = -np.log10(np.absolute(value[i]))
    #        continue
    #    else:               # 0
    #        tmp_array[i] = 0
    
    return tmp_array/1e11
    
        
def gradient_x(grid,x_value,y_value):
    return np.divide((grid[y_value,x_value]-grid[y_value,x_value-1]), (grid_x[y_value,x_value]-grid_x[y_value,x_value-1]))

def gradient_y(grid,x_value,y_value):
    return np.divide((grid[y_value,x_value]-grid[y_value-1,x_value]), (grid_y[y_value,x_value]-grid_y[y_value-1,x_value]))

gradient_gravity_x     = scale(gradient_x(grid_F_gravity,x,y))
gradient_gravity_y     = scale(gradient_y(grid_F_gravity,x,y))
gradient_gas_x         = scale(gradient_x(grid_F_gas,x,y))
gradient_gas_y         = scale(gradient_y(grid_F_gas,x,y))
gradient_magnetic_x    = scale(gradient_x(grid_F_magnetic,x,y))
gradient_magnetic_y    = scale(gradient_y(grid_F_magnetic,x,y))
gradient_centrifugal_x = gradient_x(grid_F_centrifugal,x,y)
gradient_centrifugal_y = gradient_y(grid_F_centrifugal,x,y)
gradient_centrifugal_X = scale(np.sqrt(np.square(gradient_centrifugal_x)+np.square(gradient_centrifugal_y)))
gradient_total_x       = gradient_gravity_x+gradient_gas_x+gradient_magnetic_x+gradient_centrifugal_x
gradient_total_y       = gradient_gravity_y+gradient_gas_y+gradient_magnetic_y+gradient_centrifugal_y

gradients_jet = np.array([[gradient_gravity_x[0], gradient_gravity_y[0]], 
                    [gradient_gas_x[0], gradient_gas_y[0]], 
                    [gradient_magnetic_x[0], gradient_magnetic_y[0]], 
                    [gradient_centrifugal_X[0], 0*gradient_centrifugal_y[0]],
                    [gradient_total_x[0], gradient_total_y[0]]])

gradients_wind = np.array([[gradient_gravity_x[1], gradient_gravity_y[1]], 
                    [gradient_gas_x[1], gradient_gas_y[1]], 
                    [gradient_magnetic_x[1], gradient_magnetic_y[1]], 
                    [gradient_centrifugal_X[1], 0*gradient_centrifugal_y[1]],
                    [gradient_total_x[1], gradient_total_y[1]]])
                    
gradients_disk = np.array([[gradient_gravity_x[2], gradient_gravity_y[2]], 
                    [gradient_gas_x[2], gradient_gas_y[2]], 
                    [gradient_magnetic_x[2], gradient_magnetic_y[2]], 
                    [gradient_centrifugal_X[2], 0*gradient_centrifugal_y[2]],
                    [gradient_total_x[2], gradient_total_y[2]]])


######### Save to file/Write to file ###########
datafile = '/Users/Anton/Desktop/Data/snapshots.npy'
if not os.path.exists(datafile):
    save_data = np.array([[gradients_jet, gradients_wind, gradients_disk]])
    np.save(datafile, save_data)
else:
    try:
        loaded_data = np.load(datafile)
    except IOError as error:
        print('Error in loading datafile:', error)
        
    append_data = np.array([[gradients_jet, gradients_wind, gradients_disk]])
    appended_data = np.append(loaded_data, append_data, axis=0)
    np.save(datafile, appended_data)
    

## Plot the interpolated data
plt.close('all')
plt.figure()
plt.rc('text', usetex=True)

c = plt.contour(grid_x, grid_y, np.log10(grid_rho), levels=[0] #extend='both', levels=np.linspace(-8,-4,60)
)
#plt.colorbar(c)

#plt.plot([0,1477000*150], [0,1477000*150])

q_jet = plt.quiver(1477000*X_jet, 1477000*Y_jet, gradients_jet[:,0], gradients_jet[:,1], 
    color=['blue', 'green', 'yellow', 'red', 'black'], linestyle=['solid','solid','solid','solid','dashed'], 
    width=.004, scale=100)

q_wind = plt.quiver(1477000*X_wind, 1477000*Y_wind, gradients_wind[:,0], gradients_wind[:,1], 
    color=['blue', 'green', 'yellow', 'red', 'black'], linestyle=['solid','solid','solid','solid','dashed'], 
    width=.004, scale=100)

q_disk = plt.quiver(1477000*X_disk, 1477000*Y_disk, gradients_disk[:,0], gradients_disk[:,1], 
    color=['blue', 'green', 'yellow', 'red', 'black'], linestyle=['solid','solid','solid','solid','dashed'], 
    width=.004, scale=100)

p1 = mline.Line2D([], [], color='blue', label='Gravity')
p2 = mline.Line2D([], [], color='green', label='Gas pressure')
p3 = mline.Line2D([], [], color='yellow', label='Magnetic pressure')
p4 = mline.Line2D([], [], color='red', label='Centrifugal')
p5 = mline.Line2D([], [], color='black', label='Total')

plt.legend(handles=[p1,p2,p3,p4,p5], fontsize='x-small')

plt.title('Force distribution w/o radiation pressure: ' + str(index))
plt.xlabel('$r [cm]$')
plt.ylabel('$r [cm]$')
plt.savefig(fileout)
plt.show()