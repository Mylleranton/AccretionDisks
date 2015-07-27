import matplotlib.pyplot as plt
import matplotlib.lines as mline
import matplotlib
import os
import numpy as np
from scipy.interpolate import griddata 
from christoffel_symbols import christoffel as ch
import STATIC_DTYPE as input_dtype

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
    
def dataindex(ix, iy):
    return ix + 252*iy
        
def gradient_radial(function, line):
    return (function[line+1]-function[line-1])/(r[line+1]-r[line-1])
    
    
    tmp_gradr = np.empty([len(r),1])
    for index in range(0, len(r)-1):
        if np.logical_or((r[index+1]-r[index-1]) == 0, index+1 >= len(r)):
            print('Error: Radial gradient returned 0 at index ', index)
            tmp_gradr[index] = 0
        else:
            tmp_gradr[index] = (function[index+1]-function[index-1])/(r[index+1]-r[index-1])
    return tmp_gradr
    
def gradient_theta(function, line):
    return (function[line+252]-function[line-252])/(theta[line+252]-theta[line-252])

    
    tmp_gradt = np.empty([len(theta),1])
    for index in range(0, len(theta)-1):
        if (index + 252>=len(theta)):
            print('Index out of bounds:', index)
            continue
        elif ((theta[index+252]-theta[index-252]) == 0):
            print('Error: Theta gradient returned 0 at index ',index)
            tmp_gradt[index] = 0
        else:
            tmp_gradt[index] = (function[index+252]-function[index-252])/(theta[index+252]-theta[index-252])
    return tmp_gradt
    
# Create meshgrid (in units of gravitational radii, then converted to physical)
ymin = 0
ymax = 100
xmin = 0
xmax = 100

grid_y, grid_x = np.mgrid[ymin:ymax:150j,xmin:xmax:150j]
grid_y = grid_y
grid_x = grid_x


## Loop through all simulation files
for index in range(700,701):
    if index < 1000:
        index_str = str('0' + str(index))
    else:
        index_str = str(index)
    
    #File operations
    filename = 'sim' + index_str
    filein = open('/Volumes/Seagate/4Anton/hd300a0/dt100/simavg0070-0134_rel.dat','rb')
    #filein = open('/Volumes/Seagate/4Anton/d300a0/' + filename + '.dat','rb')
    #datafile_all_forces = '/Users/Anton/Desktop/Data/Binaries/RAD-600-d300a0/' + index_str + '.npy'
    #datafile_all_forces = '/Users/Anton/Desktop/Data/Binaries/RAD-avg-d300a0/average.npy'


    #fileout = '/Users/Anton/Desktop/Data/Image/' + filename + '.png'
        
    
    datain = np.loadtxt(filein,input_dtype.dtype_rel_hydro())
    filein.close()
    #datain = datain[0:10000]
    LINE_INDEX = np.array([dataindex(200, 200), dataindex(100,20), dataindex(100,50), dataindex(100,85),  dataindex(130, 10), dataindex(130, 30),  dataindex(130, 55), dataindex(130, 85), dataindex(145, 10), dataindex(145, 30), dataindex(145, 55),  dataindex(145, 85),dataindex(155,30)])
    NUMBER_OF_FORCES = 1
   
    ########### ------------------------------------------------------------------------------------ ###########
    ########### ------------------------------------------------------------------------------------ ###########
    ########### ------------------------------------------------------------------------------------ ###########
    
    #Convert coordinates from spherical to cartesian (physical units)
    r = datain['r']
    theta = datain['theta']
    phi = datain['phi']
    rho = datain['rho']
    
    points = (r*np.sin(theta),r*np.cos(theta))
    grid_rho = griddata(points, rho, (grid_x, grid_y), method='linear')
    
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
    
    
    F_x = np.empty((NUMBER_OF_FORCES,1))
    F_y = np.empty((NUMBER_OF_FORCES,1))
    x = np.empty((NUMBER_OF_FORCES,1))
    y = np.empty((NUMBER_OF_FORCES,1))
    
    for index2 in range(0,len(LINE_INDEX)):
        line = LINE_INDEX[index2]
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
        
        corr_gradients = gradient_radial((w-rho)*u_radial, line) + gradient_theta((w-rho)*u_theta, line)
        F_enth_corr_radial = -(1./w)*u_radial*g11*corr_gradients
        F_enth_corr_theta = -(1./w)*u_theta*g22*corr_gradients/r
        
        F_correction_radial = F_rel_corr_radial + F_enth_corr_radial
        F_correction_theta = F_rel_corr_theta + F_enth_corr_theta
        
        F_correction_x = F_correction_radial*np.sin(theta) + F_correction_theta*np.cos(theta)
        F_correction_y = F_correction_radial*np.cos(theta) - F_correction_theta*np.sin(theta)
        
        #print('010', (T00*ch(r,theta,a,0,1,0)/w)[line])
        #print('111', (T11*ch(r,theta,a,1,1,1)/w)[line])
        #print('212', (T22*ch(r,theta,a,2,1,2)/w)[line])
        #print('313', (T33*ch(r,theta,a,3,1,3)/w)[line])
        #print('T11', ((1./w)*((T11-rho*u_radial*u_radial*g11)*(2./r)))[line])
        #print('T21', ((1./w)*((T21-rho*u_theta*u_radial*g11)*(np.cos(theta)/np.sin(theta))))[line])
        #print('pressure', ((2./r)*(1./w)*((5./3.-1)*datain['u_internal']+datain['bsq']/2.))[line])
        #print('enth', F_enth_corr_radial[line])
        print('Non-zero test:', (F_metric_radial[line] - (1./w[line])*((T00*ch(r,theta,a,0,1,0))[line] + (T11*ch(r,theta,a,1,1,1))[line] + (T22*ch(r,theta,a,2,1,2))[line] + (T33*ch(r,theta,a,3,1,3))[line]) + (((1./w[line])*((T11-rho*u_radial*u_radial*g11)*(2./r)))[line] + ((1./w[line])*((T21-rho*u_theta*u_radial*g11)*(np.cos(theta)/np.sin(theta))))[line])
))
  
    
        
        ## Thermal force (gas pressure) ##
        gas_pressure = ((5./3.)-1.)*datain['u_internal']
        
        F_pressure_radial = -(1./w)*gradient_radial(gas_pressure, line)
        F_pressure_theta = -(1./w)*gradient_theta(gas_pressure, line)/r
        
        F_pressure_x = (F_pressure_radial*np.sin(theta) + F_pressure_theta*np.cos(theta))
        F_pressure_y = (F_pressure_radial*np.cos(theta) - F_pressure_theta*np.sin(theta))
        
        
        
        ## Magnetic force ##
        F_magnetic_tension_radial = (1./w)*(gradient_radial(datain['b_1']*datain['b_1']*g11, line)+gradient_theta(datain['b_2']*datain['b_1']*g11, line))
        F_magnetic_tension_theta = (1./w)*(gradient_radial(datain['b_1']*datain['b_2']*g22, line)+gradient_theta(datain['b_2']*datain['b_2']*g22, line))/r
        
        F_magnetic_pressure_radial = -(1./w)*gradient_radial(datain['bsq']/2., line)
        F_magnetic_pressure_theta = -(1./w)*gradient_theta(datain['bsq']/2., line)/r

        
        F_magnetic_radial = F_magnetic_pressure_radial + F_magnetic_tension_radial
        F_magnetic_theta = F_magnetic_pressure_theta + F_magnetic_tension_theta
        
        F_magnetic_x = (F_magnetic_theta*np.cos(theta) + F_magnetic_radial*np.sin(theta))
        F_magnetic_y = (-F_magnetic_theta*np.sin(theta) + F_magnetic_radial*np.cos(theta))
        
        
        ## Total force ##
        F_total_x = F_gravity_x + F_pressure_x + F_centrifugal_x + F_magnetic_x + F_correction_x 
        F_total_y = F_gravity_y + F_pressure_y + F_centrifugal_y + F_magnetic_y + F_correction_y 
        
        if index2 == 0:
            norm_factor = np.sqrt(F_total_x[line]**(2.)+F_total_y[line]**(2.))
            F_x[0] = F_total_x[line]/norm_factor
            F_y[0] = F_total_y[line]/norm_factor
            x[0:NUMBER_OF_FORCES] = r[line]*np.sin(theta[line])
            y[0:NUMBER_OF_FORCES] = r[line]*np.cos(theta[line])
        else:
            norm_factor = np.sqrt(F_total_x[line]**(2.)+F_total_y[line]**(2.))
            tmp_Fx = np.array([[F_total_x[line]/norm_factor]] )
            tmp_Fy = np.array([[F_total_y[line]/norm_factor]] )
            
            tmp_x = (np.empty((NUMBER_OF_FORCES,1))[0:NUMBER_OF_FORCES])
            tmp_x[0:NUMBER_OF_FORCES] = r[line]*np.sin(theta[line])
            tmp_y = (np.empty((NUMBER_OF_FORCES,1))[0:NUMBER_OF_FORCES])
            tmp_y[0:NUMBER_OF_FORCES] = r[line]*np.cos(theta[line])
    
            F_x = np.append(F_x, tmp_Fx, axis=0)
            F_y = np.append(F_y, tmp_Fy, axis=0)  
            x = np.append(x, tmp_x, axis=0)
            y = np.append(y, tmp_y, axis=0) 
               
    #grid_magnetic_pressure = griddata(points, datain['bsq']/2., (grid_x, grid_y), method='linear')
    #grid_gas_pressure = griddata(points, gas_pressure, (grid_x, grid_y), method='linear')
    grid_u_magnitude = griddata(points, np.sqrt(u_radial*u_radial + u_theta*u_theta*g22*g22), (grid_x, grid_y), method='linear')
    #ux = (-u_theta*np.cos(theta)*r + u_radial*np.sin(theta))
    #uy = (-u_theta*np.sin(theta)*r + u_radial*np.cos(theta))
    #grid_ux  = griddata(points, ux, (grid_x, grid_y), method='linear')
    #grid_uy  = griddata(points, uy, (grid_x, grid_y), method='linear')
        
    fig = plt.figure()
    ax = fig.gca()
    plt.rc('text', usetex=True)
    
    ymin = 0
    ymax = 100
    xmin = 0
    xmax = 70
    
    plt.xlim(xmin=xmin, xmax=xmax)
    plt.ylim(ymin=ymin, ymax=ymax)
    
    c = plt.contourf(grid_x, grid_y, grid_u_magnitude, extend='both', levels=np.linspace(0,0.1,41), cmap='gray', alpha=0.4)
    
    cax, args = matplotlib.colorbar.make_axes(ax, location='top', shrink=0.45)
    bar = fig.colorbar(c, cax=cax, ticks=np.linspace(0,0.1,5), orientation='horizontal')
    bar.ax.xaxis.set_ticks_position('top')    
    plt.sca(ax)
    

    #x = r[LINE_INDEX]*np.sin(theta[LINE_INDEX])
    #y = r[LINE_INDEX]*np.cos(theta[LINE_INDEX])
    
    #F_x = np.array([F_gravity_x[LINE_INDEX], F_pressure_x[LINE_INDEX], F_magnetic_x[LINE_INDEX], F_centrifugal_x[LINE_INDEX], F_total_x[LINE_INDEX]])
    #F_y = np.array([F_gravity_y[LINE_INDEX], F_pressure_y[LINE_INDEX], F_magnetic_y[LINE_INDEX],                         0, F_total_y[LINE_INDEX]])
    
    plt.quiver(x,y,F_x,F_y, 
            color=['red'], scale=7.5, width=0.015)
    
    #p1 = mline.Line2D([], [], color='blue', label='Gravity')
    #p2 = mline.Line2D([], [], color='green', label='Thermal')
    #p3 = mline.Line2D([], [], color='yellow', label='Magnetic')
    #p4 = mline.Line2D([], [], color='red', label='Centrifugal')
    #p5 = mline.Line2D([], [], color='cyan', label='Rel. correction')
    #p6 = mline.Line2D([], [], color='blue', label='Net force')
    
    #plt.legend(handles=[p6])#, fontsize='normal')
    
    #lw = 1+grid_u_magnitude*30.
    #s = plt.streamplot(grid_x, grid_y, grid_ux, grid_uy, 
    #    density=0.5, color='#DDDDDD', 
    #    arrowsize=2.5, arrowstyle='->',
    #    minlength=.1, linewidth=lw)
        
    plt.gca().set_aspect('equal')    
    #plt.title('Velocity and net force')
    plt.xlabel('$x/r_g$')
    #plt.ylabel('$r/r_g$')
    plt.tick_params(axis='both', which='both', bottom='on', top='off', labelbottom='on', right='off', left='off', labelleft='off')

    fileout = '/Users/Anton/Dropbox/Aleksander/Figures/simavg0070-0134/net_force.png'
    plt.savefig(fileout, bbox_inches='tight')     
    
    plt.show()

    
    print('Index:' + index_str)