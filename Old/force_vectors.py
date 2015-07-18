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
filename = 'simavg0070-0134'
filein = open('/Users/Anton/Desktop/Data/hd300a0/' + filename + '.dat','rb')
fileout = '/Users/Anton/Desktop/Data/Image/' + filename + '_2.png'

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

        
def gradient_x(grid,x_value,y_value):
    return np.divide((grid[y_value,x_value]-grid[y_value,x_value-1]), (grid_x[y_value,x_value]-grid_x[y_value,x_value-1]))

def gradient_y(grid,x_value,y_value):
    return np.divide((grid[y_value,x_value]-grid[y_value-1,x_value]), (grid_y[y_value,x_value]-grid_y[y_value-1,x_value]))
    

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
F_gravity = 6.674*np.power(10.00, -8)*10.*1.9891*np.power(10.00, 33)/(np.square(1477000.*datain['r']))
F_gravity_x = -F_gravity*np.sin(datain['theta'])
F_gravity_y = -F_gravity*np.cos(datain['theta'])

## Thermal force (gas pressure) ##
# P = b*rho*T/(mu * protonmass)
# 1.38e-16
# 1.67e-24
# F_gas = -gradient(P)/rho
thermal_pressure = (np.divide(1.38,1.67)*np.multiply(datain['rho'],datain['temp'])*np.power(10.,8))

## Magnetic force (magnetic pressure) ##
magnetic_pressure = -np.divide(1, datain['rho'])*np.divide(1, 2)

## Centrifugal force ##
v_phi = (np.divide(datain['u_3'], datain['u_t'])*datain['r']*np.sin(datain['theta'])*3.e10)

effective_radius = 1477000.*datain['r']*np.sin(datain['theta'])
F_centrifugal = np.divide(np.square(v_phi), effective_radius)


##### Interpolate forces and vector fields #####
grid_F_gravity_x         = griddata(points, F_gravity_x, (grid_x, grid_y), method='linear')
grid_F_gravity_y        = griddata(points, F_gravity_y, (grid_x, grid_y), method='linear')

grid_thermal_pressure     = griddata(points, thermal_pressure, (grid_x, grid_y), method='linear')

grid_magnetic_pressure    = griddata(points, magnetic_pressure, (grid_x, grid_y), method='linear')
grid_bsq                = griddata(points, datain['bsq'], (grid_x, grid_y), method='linear')
grid_F_centrifugal     = griddata(points, F_centrifugal, (grid_x, grid_y), method='linear')
grid_rho               = griddata(points, datain['rho'], (grid_x, grid_y), method='linear')

########### ------------------------------------------------------------------------------------ ###########
########### ------------------------------------------------------------------------------------ ###########
########### ------------------------------------------------------------------------------------ ###########


# Wind = (40,70), Jet = (15,80), Disk = (80,20)
x         = np.array([0,40,40])
X_all     = 1477000*np.array([x[0],x[0],x[0],x[0],x[0], x[1],x[1],x[1],x[1],x[1], x[2],x[2],x[2],x[2],x[2]])
y         = np.array([0,15,2])
Y_all     = 1477000*np.array([y[0],y[0],y[0],y[0],y[0],y[1],y[1],y[1],y[1],y[1],y[2],y[2],y[2],y[2],y[2]])

#### Returns the logarithm of the value while preserving the sign.
def scale(value):
    tmp_array = value
    return tmp_array/1.e11

    

gravity_x             = scale(grid_F_gravity_x[y,x])
gravity_y             = scale(grid_F_gravity_y[y,x])

gradient_gas_x         = scale(-gradient_x(grid_thermal_pressure,x,y)/grid_rho[y,x])
gradient_gas_y         = scale(-gradient_y(grid_thermal_pressure,x,y)/grid_rho[y,x])

gradient_magnetic_x    = scale(-gradient_x(grid_bsq,x,y)/grid_rho[y,x])
gradient_magnetic_y    = scale(-gradient_y(grid_bsq,x,y)/grid_rho[y,x])

centrifugal_x         = scale(grid_F_centrifugal[y,x])
gradient_total_x       = gravity_x+gradient_gas_x+gradient_magnetic_x+centrifugal_x
gradient_total_y       = gravity_y+gradient_gas_y+gradient_magnetic_y+0

gradients_all = np.array([[gravity_x[0], gravity_y[0]], 
                    [gradient_gas_x[0], gradient_gas_y[0]], 
                    [gradient_magnetic_x[0], gradient_magnetic_y[0]], 
                    [centrifugal_x[0], 0],
                    [gradient_total_x[0], gradient_total_y[0]],
                    [gravity_x[1], gravity_y[1]], 
                    [gradient_gas_x[1], gradient_gas_y[1]], 
                    [gradient_magnetic_x[1], gradient_magnetic_y[1]], 
                    [centrifugal_x[1], 0],
                    [gradient_total_x[1], gradient_total_y[1]],
                    [gravity_x[2], gravity_y[2]], 
                    [gradient_gas_x[2], gradient_gas_y[2]], 
                    [gradient_magnetic_x[2], gradient_magnetic_y[2]], 
                    [centrifugal_x[2], 0],
                    [gradient_total_x[2], gradient_total_y[2]]])


## Plot the interpolated data

plt.figure()
plt.rc('text', usetex=True)

#ymin = 1477000*ymin
#ymax = 1477000*ymax
#xmin = 1477000*xmin
#xmax = 1477000*xmax

plt.xlim(xmin=xmin, xmax=xmax)
plt.ylim(ymin=ymin, ymax=ymax)

c = plt.contour(grid_x/1477000, grid_y/1477000, grid_rho/1477000) #, extend='both', levels=np.linspace(8,16,60))
#plt.colorbar(c)


q_all = plt.quiver(X_all/1477000, Y_all/1477000, gradients_all[:,0], gradients_all[:,1], 
    color=['blue', 'green', 'yellow', 'red', 'black'], linestyle=['solid','solid','solid','solid','dashed'], 
    width=.004, scale=20)


p1 = mline.Line2D([], [], color='blue', label='Gravity')
p2 = mline.Line2D([], [], color='green', label='Gas pressure')
p3 = mline.Line2D([], [], color='yellow', label='Magnetic pressure')
p4 = mline.Line2D([], [], color='red', label='Centrifugal')
p5 = mline.Line2D([], [], color='black', label='Total')

plt.legend(handles=[p1,p2,p3,p4,p5], fontsize='x-small')

plt.annotate(s='Jet (' + str(x[0]) + ',' + str(y[0]) + ')', xy=(x[0],y[0]), 
    size='x-small', xytext=(50, 30), textcoords='offset points')
plt.annotate(s='Wind (' + str(x[1]) + ',' + str(y[1]) + ')', xy=(x[1],y[1]),
    size='x-small', xytext=(0, -50), textcoords='offset points')
plt.annotate(s='Disk (' + str(x[2]) + ',' + str(y[2]) + ')', xy=(x[2],y[2]),
    size='x-small', xytext=(0, 30), textcoords='offset points')

plt.title('Force distribution and density')
plt.xlabel('$r/r_g$')
plt.ylabel('$r/r_g$')
plt.savefig(fileout)
plt.show()