import matplotlib.pyplot as plt
import matplotlib.lines as mline
import matplotlib
import os
import numpy as np
from christoffel_symbols import christoffel as ch
import force_vectors_all_forces as fv 

matplotlib.rcdefaults()
matplotlib.rc('font',**{'family':'serif','size':16})
os.environ['PATH'] = os.environ['PATH'] + ':/opt/local/bin'
os.environ['PATH']



def dataindex(x_loc,y_loc,datain_loc):
    x_loc = float(x_loc)
    y_loc = float(y_loc)
    radius = np.sqrt(x_loc**(2.)+y_loc**(2.))
    angle = np.arctan(x_loc/y_loc)
    line_radius = datain_loc[np.nanargmin((np.abs(datain_loc['r']-(radius))), axis=0)]
    radius = line_radius[3]
    line_angle = datain_loc[datain_loc['r']==radius][np.nanargmin((np.abs((datain_loc[datain_loc['r']==radius]['theta'])-(angle))), axis=0)]
    angle = line_angle[4]
    return np.nanargmin((np.abs(datain_loc['r']-(radius)) + np.abs(datain_loc ['theta']-(angle))), axis=0)
    
def plotforces():
     plt.quiver(x,y,F_x,F_y, 
            color=['blue', 'green', 'darkmagenta', 'red', 'cyan', 'black'], scale=.005)
     plt.show()

def plotarrows():
    plt.figure()
    plt.rc('text', usetex=True)
    ypmin = 0
    ypmax = 100
    xpmin = 0
    xpmax = 70
        
    plt.xlim(xmin=xpmin, xmax=xpmax)
    plt.ylim(ymin=ypmin, ymax=ypmax)
    
    for i in range(0,len(COORDINATES)-1):
        plt.arrow(COORDINATES[i,0], COORDINATES[i,1], (COORDINATES[i+1,0] - COORDINATES[i,0]), (COORDINATES[i+1,1] - COORDINATES[i,1]))
        #if i%8==0:
        #    plt.arrow(COORDINATES[i,0], COORDINATES[i,1], (COORDINATES[i+1,0] - COORDINATES[i,0]), (COORDINATES[i+1,1] - COORDINATES[i,1]), fc="k", ec="k", head_width=1.5, head_length=1)
        #else:
        #    plt.arrow(COORDINATES[i,0], COORDINATES[i,1], (COORDINATES[i+1,0] - COORDINATES[i,0]), (COORDINATES[i+1,1] - COORDINATES[i,1]))
    
    plt.gca().set_aspect('equal')    
    plt.title('Particle trail')
    plt.xlabel('$r/r_g$')
    plt.ylabel('$r/r_g$')


ymin = 0.#*1477000.
ymax = 100.#*1477000.
xmin = 0.#*1477000.
xmax = 100.#*1477000.
resolution = 100.
savefilename = '/Users/Anton/Desktop/Data/Binaries/hydro_particle_7_50.npy'
    
    
#File operations
filebase = '/Users/Anton/Desktop/Data/hd300a0/hd300a0_rel/'
filenumber_start = 1000
filenumber = filenumber_start

try:
    filein = open(filebase + 'sim' + str(filenumber) + '.dat','rb')
    #filein = open(filebase + 'simavg0070-0134_rel.dat', 'rb')
    datain = np.loadtxt(filein,fv.input_dtype.dtype_rel_hydro())
    filein.close()
except IndexError:
    print('Wrong input dtype. Terminating.')
    exit()
    
COORDINATES = np.load(savefilename)
step = len(COORDINATES)/4.
coord_index = np.floor(np.array([0, step, 2.*step, 3.*step, 4*step-1]))
LINE_INDEX = np.empty((len(coord_index), 1))
for i in range(0,len(coord_index)):
    tmp = COORDINATES[coord_index[i]]
    LINE_INDEX[i] = dataindex(tmp[0],tmp[1], datain)
    
    
NUMBER_OF_FORCES = 6    
F_x = np.empty((NUMBER_OF_FORCES,1))
F_y = np.empty((NUMBER_OF_FORCES,1))
x = np.empty((NUMBER_OF_FORCES,1))
y = np.empty((NUMBER_OF_FORCES,1))


for index3 in range(0,len(coord_index)):
    filenumber = filenumber_start + int(coord_index[index3])
    filein = open(filebase + 'sim' + str(filenumber) + '.dat','rb')
    #print('Loaded file: ', filenumber, 'with index: ', index3)
    #continue
    datain = np.loadtxt(filein,fv.input_dtype.dtype_rel_hydro())
    filein.close()   

    ########### ------------------------------------------------------------------------------------ ###########
        
    #Convert coordinates from spherical to cartesian (physical units)
    r = datain['r']
    theta = datain['theta']
    phi = datain['phi']
    rho = datain['rho']
        
    points = (r*np.sin(theta),r*np.cos(theta))
    #grid_rho = griddata(points, rho, (grid_x, grid_y), method='linear')
        
    ## Metric determinant ##
    a=0
    g_sch = np.power(r, 2)*np.sin(theta)
        
    g00 = -(1-np.divide(2.*r,(np.power(r, 2)+a**(2.)*np.power(np.cos(theta),2))))
    g11 = np.sqrt((np.divide( (np.power(r, 2)+a**(2.)*np.power(np.cos(theta),2)), (np.power(r, 2)-2.*r+a**2) ) ))
    g22 = np.sqrt((np.power(r,2)+a**(2.)*np.power(np.cos(theta),2)))
    g33 = np.sqrt(np.power(np.sin(theta), 2)*(   np.power(r, 2)  + a**(2.) + np.divide(2.*r, (np.power(r, 2)-2.*r+a**(2.)))*a**(2.)*np.power(np.sin(theta),2)))
        
    ##
    T00 = datain['T_00']
    T01 = datain['T_01']
    T02 = datain['T_02']
    T03 = datain['T_03']
    T10 = datain['T_10']
    T11 = datain['T_11']
    T12 = datain['T_12']
    T13 = datain['T_13']
    T20 = datain['T_20']
    T21 = datain['T_21']
    T22 = datain['T_22']
    T23 = datain['T_23']
    T30 = datain['T_30']
    T31 = datain['T_31']
    T32 = datain['T_32']
    T33 = datain['T_33']
        
        
    ## Enthalpy ##
    w = rho+datain['u_internal']+((5./3.)-1)*datain['u_internal']+datain['bsq']
    
    ## Velocities ##
    u_radial = (datain['u_1']/datain['u_t'])
    u_theta = (datain['u_2']/datain['u_t'])
    u_phi = (datain['u_3']/datain['u_t'])
        
    line = int(LINE_INDEX[index3,0])
    
    ####### Forces ########
    ## Gravitational force ##
    F_metric_radial     = (1./w)*(T00*ch(r,theta,a,0,1,0) + T01*ch(r,theta,a,1,1,0) + T02*ch(r,theta,a,2,1,0) + T03*ch(r,theta,a,3,1,0) +
                                    T10*ch(r,theta,a,0,1,1) + T11*ch(r,theta,a,1,1,1) + T12*ch(r,theta,a,2,1,1) + T13*ch(r,theta,a,3,1,1) +
                                    T20*ch(r,theta,a,0,1,2) + T21*ch(r,theta,a,1,1,2) + T22*ch(r,theta,a,2,1,2) + T23*ch(r,theta,a,3,1,2) +
                                    T30*ch(r,theta,a,0,1,3) + T31*ch(r,theta,a,1,1,3) + T32*ch(r,theta,a,2,1,3) + T33*ch(r,theta,a,3,1,3)) - (
                                    (1./w)*((T11-rho*u_radial*u_radial*g11)*(2./r)) + (1./w)*((T21-rho*u_theta*u_radial*g11)*(np.cos(theta)/np.sin(theta))))
    
    F_metric_theta     = ((1./w)*( T00*ch(r,theta,a,0,2,0) + T01*ch(r,theta,a,1,2,0) + T02*ch(r,theta,a,2,2,0) + T03*ch(r,theta,a,3,2,0) +
                                    T10*ch(r,theta,a,0,2,1) + T11*ch(r,theta,a,1,2,1) + T12*ch(r,theta,a,2,2,1) + T13*ch(r,theta,a,3,2,1) +
                                    T20*ch(r,theta,a,0,2,2) + T21*ch(r,theta,a,1,2,2) + T22*ch(r,theta,a,2,2,2) + T23*ch(r,theta,a,3,2,2) +
                                    T30*ch(r,theta,a,0,2,3) + T31*ch(r,theta,a,1,2,3) + T32*ch(r,theta,a,2,2,3) + T33*ch(r,theta,a,3,2,3)) - (
                                    (1./w)*((T12-rho*u_radial*u_theta*g22)*(2./r)) + (1./w)*((T22-rho*u_theta*u_theta*g22)*(np.cos(theta)/np.sin(theta)))))/r
                                    
    F_gravity_radial     = (1./w)*(T00-0*((5./3.)-1.)*datain['u_internal']-datain['bsq']/2.)*ch(r,theta,a,0,1,0)
    F_gravity_theta      = 0 
    F_gravity_x = F_gravity_radial*np.sin(theta)
    F_gravity_y = F_gravity_radial*np.cos(theta) 
    F_centrifugal_radial = (1./w)*T33*ch(r,theta,a,3,1,3)
    F_centrifugal_theta = (1./w)*T33*ch(r,theta,a,3,2,3)/r 
    F_centrifugal_x = F_centrifugal_radial*np.sin(theta) + F_centrifugal_theta*np.cos(theta)
    F_centrifugal_y = F_centrifugal_radial*np.cos(theta) - F_centrifugal_theta*np.sin(theta) 
        
        
    ## Relativistic correction
    F_rel_corr_radial = F_metric_radial - F_gravity_radial - F_centrifugal_radial 
    F_rel_corr_theta = F_metric_theta - F_gravity_theta - F_centrifugal_theta  
    corr_gradients = fv.gradient_radial((w-rho)*u_radial, line) + fv.gradient_theta((w-rho)*u_theta, line)
    F_enth_corr_radial = -(1./w)*u_radial*g11*corr_gradients
    F_enth_corr_theta = -(1./w)*u_theta*g22*corr_gradients/r
    F_correction_radial = F_rel_corr_radial + F_enth_corr_radial
    F_correction_theta = F_rel_corr_theta + F_enth_corr_theta
    F_correction_x = F_correction_radial*np.sin(theta) + F_correction_theta*np.cos(theta)
    F_correction_y = F_correction_radial*np.cos(theta) - F_correction_theta*np.sin(theta)
    print('Non-zero test:', (F_metric_radial[line] - (1./w[line])*((T00*ch(r,theta,a,0,1,0))[line] + (T11*ch(r,theta,a,1,1,1))[line] + (T22*ch(r,theta,a,2,1,2))[line] + (T33*ch(r,theta,a,3,1,3))[line]) + (((1./w[line])*((T11-rho*u_radial*u_radial*g11)*(2./r)))[line] + ((1./w[line])*((T21-rho*u_theta*u_radial*g11)*(np.cos(theta)/np.sin(theta))))[line])
    ))
    ## Thermal force (gas pressure) ##
    gas_pressure = ((5./3.)-1.)*datain['u_internal']
    F_pressure_radial = -(1./w)*fv.gradient_radial(gas_pressure, line)
    F_pressure_theta = -(1./w)*fv.gradient_theta(gas_pressure, line)/r
    F_pressure_x = (F_pressure_radial*np.sin(theta) + F_pressure_theta*np.cos(theta))
    F_pressure_y = (F_pressure_radial*np.cos(theta) - F_pressure_theta*np.sin(theta))
    
    ## Magnetic force ##
    F_magnetic_tension_radial = (1./w)*(fv.gradient_radial(datain['b_1']*datain['b_1']*g11, line)+fv.gradient_theta(datain['b_2']*datain['b_1']*g11, line))
    F_magnetic_tension_theta = (1./w)*(fv.gradient_radial(datain['b_1']*datain['b_2']*g22, line)+fv.gradient_theta(datain['b_2']*datain['b_2']*g22, line))/r
    F_magnetic_pressure_radial = -(1./w)*fv.gradient_radial(datain['bsq']/2., line)
    F_magnetic_pressure_theta = -(1./w)*fv.gradient_theta(datain['bsq']/2., line)/r    
    F_magnetic_radial = F_magnetic_pressure_radial + F_magnetic_tension_radial
    F_magnetic_theta = F_magnetic_pressure_theta + F_magnetic_tension_theta
    F_magnetic_x = (F_magnetic_theta*np.cos(theta) + F_magnetic_radial*np.sin(theta))
    F_magnetic_y = (-F_magnetic_theta*np.sin(theta) + F_magnetic_radial*np.cos(theta))
    
    ## Total force ##
    F_total_x = F_gravity_x + F_pressure_x + F_centrifugal_x + F_magnetic_x + F_correction_x 
    F_total_y = F_gravity_y + F_pressure_y + F_centrifugal_y + F_magnetic_y + F_correction_y
    
    if index3 == 0:
        F_x[0] = F_gravity_x[line]
        F_x[1] = F_pressure_x[line]
        F_x[2] = F_magnetic_x[line]
        F_x[3] = F_centrifugal_x[line]
        F_x[4] = F_correction_x[line]
        F_x[NUMBER_OF_FORCES-1] = F_total_x[line]
        F_y[0] = F_gravity_y[line]
        F_y[1] = F_pressure_y[line]
        F_y[2] = F_magnetic_y[line]
        F_y[3] = F_centrifugal_y[line]
        F_y[4] = F_correction_y[line]
        F_y[NUMBER_OF_FORCES-1] = F_total_y[line]
        print('LOOP1: F_x', np.shape(F_x), F_x, 'Index: ', index3)
        
        x[0:NUMBER_OF_FORCES] = r[line]*np.sin(theta[line])
        y[0:NUMBER_OF_FORCES] = r[line]*np.cos(theta[line])
    else:
        tmp_Fx = np.array([[F_gravity_x[line]], [F_pressure_x[line]], [F_magnetic_x[line]], [F_centrifugal_x[line]], [F_correction_x[line]], [F_total_x[line]]] )
        tmp_Fy = np.array([[F_gravity_y[line]], [F_pressure_y[line]], [F_magnetic_y[line]], [F_centrifugal_y[line]], [F_correction_y[line]], [F_total_y[line]]] )
        
        tmp_x = (np.empty((NUMBER_OF_FORCES,1))[0:NUMBER_OF_FORCES])
        tmp_x[0:NUMBER_OF_FORCES] = r[line]*np.sin(theta[line])
        tmp_y = (np.empty((NUMBER_OF_FORCES,1))[0:NUMBER_OF_FORCES])
        tmp_y[0:NUMBER_OF_FORCES] = r[line]*np.cos(theta[line])
    
        F_x = np.append(F_x, tmp_Fx, axis=0)
        F_y = np.append(F_y, tmp_Fy, axis=0)  
        x = np.append(x, tmp_x, axis=0)
        y = np.append(y, tmp_y, axis=0) 
        print('LOOP2 F_x', np.shape(F_x), F_x, 'Index: ', index3)
        
        
plotarrows()
plotforces()
        