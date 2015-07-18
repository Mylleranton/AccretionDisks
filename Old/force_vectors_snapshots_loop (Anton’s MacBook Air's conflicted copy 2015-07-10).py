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
                       
def scale(value):
        tmp_array = value
        return tmp_array/1.e11
    
## Returns the gradient in the x and y direction respectivly for a given grid at 
## a given coordinate
def gradient_x(grid, x_array, y_array):
    tmp_gradx = np.zeros((len(y_array),len(x_array)))
    for i in range(0, len(x_array)):
        for k in range(0,len(y_array)):
            tmp_gradx[k,i] = np.divide(np.subtract(grid[y_array[k], x_array[i]], grid[y_array[k], x_array[i]-1]), np.subtract(grid_x[y_array[k], x_array[i]], grid_x[y_array[k], x_array[i]-1]))
            
    return tmp_gradx

def gradient_y(grid,x_array,y_array):
    tmp_grady = np.zeros((len(y_array),len(x_array)))
    for i in range(0, len(x_array)):
        for k in range(0,len(y_array)):
            tmp_grady[k,i] = np.divide(np.subtract(grid[y_array[k], x_array[i]], grid[y_array[k]-1, x_array[i]]), np.subtract(grid_y[y_array[k], x_array[i]], grid_y[y_array[k]-1, x_array[i]]))
            
    return tmp_grady
   # return np.divide((grid[y_value,x_value]-grid[y_value-1,x_value]), (grid_y[y_value,x_value]-grid_y[y_value-1,x_value]))

# Create meshgrid (in units of gravitational radii, then converted to physical)
ymin = 0
ymax = 100
xmin = 0
xmax = 100

grid_y, grid_x = np.mgrid[ymin:ymax:101j,xmin:xmax:101j]
grid_y = 1477000.*grid_y
grid_x = 1477000.*grid_x

## Plotting utilities
# Wind = (40,70), Jet = (15,80), Disk = (80,20)
x         = np.array([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100])
y         = np.array([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100])
X_all     = np.array([[x[0],x[0],x[0],x[0],x[0]]], dtype=np.dtype(int))
Y_all     = np.array([[y[0],y[0],y[0],y[0],y[0]]], dtype=np.dtype(int))

for i in range(1,len(x)):
    X_all = np.append(X_all, np.array([[x[i],x[i],x[i],x[i],x[i]]]), axis=0)

for i in range(1,len(y)):
    Y_all = np.append(Y_all, np.array([[y[i],y[i],y[i],y[i],y[i]]]), axis=0)

X_all = 1477000.*np.ndarray.flatten(X_all)
Y_all = 1477000.*np.ndarray.flatten(Y_all)


## Loop through all simulation files
for index in range(700,701):
    if index < 1000:
        index_str = str('0' + str(index))
    else:
        index_str = str(index)
    
    #File operations
    #filename = 'sim' + index_str
    filename = 'simavg0070-0134_rel'
    filein = open('/Users/Anton/Desktop/Data/hd300a0/' + filename + '.dat','rb')
    #filein = open('/Volumes/Seagate/4Anton/hd300a0/' + filename + '.dat','rb')
    #datafile_all_forces = '/Users/Anton/Desktop/Data/Binaries/NO-RAD-600-hd300a0/' + index_str + '.npy'
    #datafile_all_forces = '/Users/Anton/Desktop/Data/Binaries/NO-RAD-600-hd300a0/average.npy'

    #fileout = '/Users/Anton/Desktop/Data/Image/' + filename + '.png'
    
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
    
    input_dtype_rel_hydro = np.dtype([
        ("ix", int),
        ("iy", int),
        ("iz", int),
        ("r", float),
        ("theta", float),
        ("phi", float),
        ("rho", float),
        ("u_internal", float),
        ("u_t", float),
        ("u_1", float),
        ("u_2", float),
        ("u_3", float),
        ("volume", float),
        ("bsq", float),
        ("b_1", float),
        ("b_2", float),
        ("b_3", float),
        ("T_00", float),
        ("T_01", float),
        ("T_02", float),
        ("T_03", float),
        ("T_10", float),
        ("T_11", float),
        ("T_12", float),
        ("T_13", float),
        ("T_20", float),
        ("T_21", float),
        ("T_22", float),
        ("T_23", float),
        ("T_30", float),
        ("T_31", float),
        ("T_32", float),
        ("T_33", float)
        ])
    
    datain = np.loadtxt(filein,input_dtype_rel_hydro)
    filein.close()
    
    ########### ------------------------------------------------------------------------------------ ###########
    ########### ------------------------------------------------------------------------------------ ###########
    ########### ------------------------------------------------------------------------------------ ###########
    
    #Convert coordinates from spherical to cartesian (physical units)
    rho = datain['rho']*6.173e15
    temp = (datain['u_internal']/datain['rho'])*7.259e12
    points = (1477000.*datain['r']*np.sin(datain['theta']),1477000.*datain['r']*np.cos(datain['theta']))
    grid_rho = griddata(points, rho, (grid_x, grid_y), method='linear')
    
    ####### Forces ########
    ## Gravitational force ##
    #
    # G = 6.674e-8 cm3/g/s2
    # M-sun = 1.9891e33 g
    F_gravity               = 6.674*np.power(10.00, -8)*10.*1.9891*np.power(10.00, 33)/(np.square(1477000.*datain['r']))
    F_gravity_x             = -F_gravity*np.sin(datain['theta'])
    F_gravity_y             = -F_gravity*np.cos(datain['theta'])
    grid_F_gravity_x        = griddata(points, F_gravity_x, (grid_x, grid_y), method='linear')
    grid_F_gravity_y        = griddata(points, F_gravity_y, (grid_x, grid_y), method='linear')
    
    gravity_x             = scale((grid_F_gravity_x[y,:])[:,x])
    gravity_y             = scale((grid_F_gravity_y[y,:])[:,x])
    
    ## Thermal force (gas pressure) ##
    # P = b*rho*T/(mu * protonmass)
    # 1.38e-16
    # 1.67e-24
    # F_gas = -gradient(P)/rho
    
    thermal_pressure         = (np.divide(1.38,1.67)*np.multiply(rho,temp)*np.power(10.,8))
    grid_thermal_pressure    = griddata(points, thermal_pressure, (grid_x, grid_y), method='linear')
    
    gradient_gas_x         = scale(-gradient_x(grid_thermal_pressure,x,y)/grid_rho[y,:][:,x])
    gradient_gas_y         = scale(-gradient_y(grid_thermal_pressure,x,y)/grid_rho[y,:][:,x])
    
    ## Magnetic force (magnetic pressure) ##
    magnetic_pressure          = -np.divide(1, rho)*np.divide(1, 2)
    grid_magnetic_pressure     = griddata(points, magnetic_pressure, (grid_x, grid_y), method='linear')
    grid_bsq                   = griddata(points, datain['bsq'], (grid_x, grid_y), method='linear')
    
    gradient_magnetic_x    = scale(-gradient_x(grid_bsq,x,y)/grid_rho[y,:][:,x])
    gradient_magnetic_y    = scale(-gradient_y(grid_bsq,x,y)/grid_rho[y,:][:,x])

    
    ## Centrifugal force ##
    v_phi                 = (np.divide(datain['u_3'], datain['u_t'])*datain['r']*np.sin(datain['theta'])*3.e10)
    effective_radius      = 1477000.*datain['r']*np.sin(datain['theta'])
    F_centrifugal         = np.divide(np.square(v_phi), effective_radius)
    grid_F_centrifugal    = griddata(points, F_centrifugal, (grid_x, grid_y), method='linear')
    centrifugal_x         = scale(grid_F_centrifugal[y,:][:,x])
   
    
    ## Total force ##
    gradient_total_x       = gravity_x+gradient_gas_x+gradient_magnetic_x+centrifugal_x
    gradient_total_y       = gravity_y+gradient_gas_y+gradient_magnetic_y+0
    
                        
    output_array = np.array([[x[0], y[0], gravity_x[0,0], gravity_y[0,0], gradient_gas_x[0,0], gradient_gas_y[0,0], gradient_magnetic_x[0,0], gradient_magnetic_y[0,0], centrifugal_x[0,0], 0, gradient_total_x[0,0], gradient_total_y[0,0]]])
    for i in range(0,len(x)):
        for k in range(1,len(y)):
            output_array = np.append(output_array, np.array([[x[i], y[k], gravity_x[k,i], gravity_y[k,i], gradient_gas_x[k,i], gradient_gas_y[k,i], gradient_magnetic_x[k,i], gradient_magnetic_y[k,i], centrifugal_x[k,i], 0, gradient_total_x[k,i], gradient_total_y[k,i]]]), axis=0)
	
        
    ######### Save to file/Write to file ###########
    #if not os.path.exists(datafile_all_forces):
    #    np.save(datafile_all_forces, output_array)
    #else:
    #    try:
    #        loaded_data = np.load(datafile_all_forces)
    #    except IOError as error:
    #        print('Error in loading datafile:', error)
    #        break
    #    np.save(datafile_all_forces, output_array)
    
    
    plt.figure()
    plt.rc('text', usetex=True)

    ymin = 0
    ymax = 100
    xmin = 0
    xmax = 100
    
    plt.xlim(xmin=xmin, xmax=xmax)
    plt.ylim(ymin=ymin, ymax=ymax)
    
    LINE_INDEX = 25320
    x = datain['r'][LINE_INDEX]*np.sin(datain['theta'][LINE_INDEX])
    y = datain['r'][LINE_INDEX]*np.cos(datain['theta'][LINE_INDEX])
    
    F_x = np.array([gravity_x[LINE_INDEX], gradient_gas_x[LINE_INDEX]])
    F_y = np.array([gravity_y[LINE_INDEX], gradient_gas_y[LINE_INDEX]])
    
    
    q_all = plt.quiver(X_all, Y_all, F_x, F_y, 
        color=['blue', 'green'])
    plt.show()
    
        
    print('Index:' + index_str)