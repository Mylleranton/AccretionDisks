import matplotlib.pyplot as plt
import matplotlib.lines as mline
import matplotlib
import os
import numpy as np
from scipy.interpolate import griddata 
from christoffel_symbols import christoffel as ch

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

grid_y, grid_x = np.mgrid[ymin:ymax:101j,xmin:xmax:101j]
grid_y = grid_y
grid_x = grid_x

## Plotting utilities
# Wind = (40,70), Jet = (15,80), Disk = (80,20)
#x         = np.array([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100])
#y         = np.array([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100])
#X_all     = np.array([[x[0],x[0],x[0],x[0],x[0],x[0]]], dtype=np.dtype(int))
#Y_all     = np.array([[y[0],y[0],y[0],y[0],y[0],y[0]]], dtype=np.dtype(int))
#
#for i in range(1,len(x)):
#    X_all = np.append(X_all, np.array([[x[i],x[i],x[i],x[i],x[i],x[i]]]), axis=0)
#
#for i in range(1,len(y)):
#    Y_all = np.append(Y_all, np.array([[y[i],y[i],y[i],y[i],y[i],y[i]]]), axis=0)
#
#X_all = 1477000.*np.ndarray.flatten(X_all)
#Y_all = 1477000.*np.ndarray.flatten(Y_all)


## Loop through all simulation files
for index in range(700,701):
    if index < 1000:
        index_str = str('0' + str(index))
    else:
        index_str = str(index)
    
    #File operations
    filename = 'sim' + index_str
    filein = open('/Users/Anton/Desktop/Data/hd300a0/simavg0070-0134_rel.dat','rb')
    #filein = open('/Volumes/Seagate/4Anton/d300a0/' + filename + '.dat','rb')
    #datafile_all_forces = '/Users/Anton/Desktop/Data/Binaries/RAD-600-d300a0/' + index_str + '.npy'
    #datafile_all_forces = '/Users/Anton/Desktop/Data/Binaries/RAD-avg-d300a0/average.npy'


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
        ("rad_ehat", float),
        ("radflux_1", float),
        ("radflux_2", float),
        ("radflux_3", float),
        ("rad_tot_en", float),
        ("rad_compt_en", float),
        ("rad_bb_temp", float),
        ("rad_bb_temp2", float),
        ("bsq", float),
        ("b_1", float),
        ("b_2", float),
        ("b_3", float)
        ])
    
    input_dtype_rel_rad = np.dtype([
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
        ("T_33", float),
        ("rad_ehat", float),
        ("R_00", float),
        ("R_01", float),
        ("R_02", float),
        ("R_03", float),
        ("R_10", float),
        ("R_11", float),
        ("R_12", float),
        ("R_13", float),
        ("R_20", float),
        ("R_21", float),
        ("R_22", float),
        ("R_23", float),
        ("R_30", float),
        ("R_31", float),
        ("R_32", float),
        ("R_33", float),
        ("G_time", float),
        ("G_radial", float),
        ("G_theta", float),
        ("G_phi", float)
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
    #datain = datain[0:10000]
    LINE_INDEX = np.array([dataindex(100,10),  dataindex(100,30), dataindex(100,55), dataindex(100,85),  dataindex(130, 10), dataindex(130, 30),  dataindex(130, 55), dataindex(130, 85), dataindex(150, 10), dataindex(150, 30), dataindex(150, 55),  dataindex(150, 85)])
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
    
    
    F_x = np.empty((6,1))
    F_y = np.empty((6,1))
    x = np.empty((6,1))
    y = np.empty((6,1))
    
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
        #print('Non-zero test:', (F_metric_radial-T00*ch(r,theta,a,0,1,0)/w-T11*ch(r,theta,a,1,1,1)/w-T22*ch(r,theta,a,2,1,2)/w-T33*ch(r,theta,a,3,1,3)/w+(1./w)*(T11-rho*u_radial*u_radial*g11*(2./r))+(1./w)*(T21-rho*u_theta*u_radial*g11*(np.cos(theta)/np.sin(theta))))[line])
        #
        #    
    
        
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
        
        
        ## Radiative force ##
        
        #F_radiation_radial = (1./w)*datain['G_radial']
        #F_radiation_theta = (1./w)*datain['G_theta']
        #
        #
        
        
        
        
        #radflux_x = -(datain['radflux_2']*np.cos(theta)*np.sqrt(g22) + datain['radflux_1']*np.sin(theta)*np.sqrt(g11))
        #radflux_y = -(-datain['radflux_2']*np.sin(theta)*np.sqrt(g22) + datain['radflux_1']*np.cos(theta)*np.sqrt(g11))
        #
        #F_radflux_x = datain['u_t']*0.34*radflux_x/3.e10
        #F_radflux_y = datain['u_t']*0.34*radflux_y/3.e10
        #
        #F_raddrag_x = -datain['u_t']*0.34*((4/3)*datain['rad_ehat']*ux)/3.e10
        #F_raddrag_y = -datain['u_t']*0.34*((4/3)*datain['rad_ehat']*ux)/3.e10
        #
        #grid_F_rad_x = griddata(points, (F_radflux_x+F_raddrag_x), (grid_x, grid_y), method='linear')
        #grid_F_rad_y = griddata(points, (F_radflux_y+F_raddrag_x), (grid_x, grid_y), method='linear')  
        #
        #radiation_x = scale(grid_F_rad_x[y,:][:,x])
        #radiation_y = scale(grid_F_rad_y[y,:][:,x])
        
        ## Total force ##
        F_total_x = F_gravity_x + F_pressure_x + F_centrifugal_x + F_magnetic_x + F_correction_x #+ F_correction_radial + F_radiation_radial
        F_total_y = F_gravity_y + F_pressure_y + F_centrifugal_y + F_magnetic_y + F_correction_y #+ F_correction_theta + F_radiation_theta
        
        if index2 == 0:
            F_x[0] = F_gravity_x[line]
            F_x[1] = F_pressure_x[line]
            F_x[2] = F_magnetic_x[line]
            F_x[3] = F_centrifugal_x[line]
            F_x[4] = F_correction_x[line]
            F_x[5] = F_total_x[line]
            F_y[0] = F_gravity_y[line]
            F_y[1] = F_pressure_y[line]
            F_y[2] = F_magnetic_y[line]
            F_y[3] = F_centrifugal_y[line]
            F_x[4] = F_correction_y[line]
            F_y[5] = F_total_y[line]
            x[0:6] = r[line]*np.sin(theta[line])
            y[0:6] = r[line]*np.cos(theta[line])
        else:
            tmp_Fx = np.array([[F_gravity_x[line]], [F_pressure_x[line]], [F_magnetic_x[line]], [F_centrifugal_x[line]], [F_correction_x[line]], [F_total_x[line]]] )
            tmp_Fy = np.array([[F_gravity_y[line]], [F_pressure_y[line]], [F_magnetic_y[line]], [F_centrifugal_y[line]], [F_correction_y[line]], [F_total_y[line]]] )
            
            tmp_x = (np.empty((6,1))[0:6])
            tmp_x[0:6] = r[line]*np.sin(theta[line])
            tmp_y = (np.empty((6,1))[0:6])
            tmp_y[0:6] = r[line]*np.cos(theta[line])
    
            F_x = np.append(F_x, tmp_Fx, axis=0)
            F_y = np.append(F_y, tmp_Fy, axis=0)  
            x = np.append(x, tmp_x, axis=0)
            y = np.append(y, tmp_y, axis=0) 
        
                    
 #   output_array = np.zeros((1,14))
 #   for i in range(0,len(x)):
 #       for k in range(0,len(y)):
 #           if (np.logical_and(i==0, k==0)):
 #               output_array = np.array([[x[0], y[0], gravity_x[0,0], gravity_y[0,0], gradient_gas_x[0,0], gradient_gas_y[0,0], F_magnetic_x[0,0], F_magnetic_y[0,0], centrifugal_x[0,0], 0, radiation_x[0,0], radiation_y[0,0],gradient_total_x[0,0], gradient_total_y[0,0]]])
 #               continue
 #           else:
 #               output_array = np.append(output_array, np.array([[x[i], y[k], gravity_x[k,i], gravity_y[k,i], gradient_gas_x[k,i], gradient_gas_y[k,i], F_magnetic_x[k,i], F_magnetic_y[k,i], centrifugal_x[k,i], 0, radiation_x[k,i], radiation_y[k,i], gradient_total_x[k,i], gradient_total_y[k,i]]]), axis=0)
	#
 #       
 #   ######### Save to file/Write to file ###########
 #   if not os.path.exists(datafile_all_forces):
 #       np.save(datafile_all_forces, output_array)
 #   else:
 #       try:
 #           loaded_data = np.load(datafile_all_forces)
 #       except IOError as error:
 #           print('Error in loading datafile:', error)
 #           break
 #       np.save(datafile_all_forces, output_array)
 #       
    #grid_magnetic_pressure = griddata(points, datain['bsq']/2., (grid_x, grid_y), method='linear')
    #grid_gas_pressure = griddata(points, gas_pressure, (grid_x, grid_y), method='linear')
    grid_u_magnitude = griddata(points, np.sqrt(u_radial*u_radial + u_theta*u_theta*g22*g22), (grid_x, grid_y), method='linear')
    
    plt.figure()
    plt.rc('text', usetex=True)
    
    ymin = 0
    ymax = 100
    xmin = 0
    xmax = 70
    
    plt.xlim(xmin=xmin, xmax=xmax)
    plt.ylim(ymin=ymin, ymax=ymax)
    
    c = plt.contourf(grid_x, grid_y, grid_u_magnitude, extend='both', levels=np.linspace(0,0.1,20), cmap='YlOrBr_r')
    plt.colorbar(c)
    

    #x = r[LINE_INDEX]*np.sin(theta[LINE_INDEX])
    #y = r[LINE_INDEX]*np.cos(theta[LINE_INDEX])
    
    #F_x = np.array([F_gravity_x[LINE_INDEX], F_pressure_x[LINE_INDEX], F_magnetic_x[LINE_INDEX], F_centrifugal_x[LINE_INDEX], F_total_x[LINE_INDEX]])
    #F_y = np.array([F_gravity_y[LINE_INDEX], F_pressure_y[LINE_INDEX], F_magnetic_y[LINE_INDEX],                         0, F_total_y[LINE_INDEX]])
    
    plt.quiver(x,y,F_x,F_y, 
            color=['blue', 'green', 'yellow', 'red', 'cyan', 'black'], scale=.005)
    
    p1 = mline.Line2D([], [], color='blue', label='Gravity')
    p2 = mline.Line2D([], [], color='green', label='Thermal')
    p3 = mline.Line2D([], [], color='yellow', label='Magnetic')
    p4 = mline.Line2D([], [], color='red', label='Centrifugal')
    p5 = mline.Line2D([], [], color='cyan', label='Rel. correction')
    p6 = mline.Line2D([], [], color='black', label='Total')
    
    plt.legend(handles=[p1,p2,p3,p4,p5,p6], fontsize='x-small')
    
      
    #s = plt.streamplot(grid_x, grid_y, grid_radflux_x, grid_radflux_y, 
    #    density=1, color='#BBBBBB', 
    #    arrowsize=2.5, arrowstyle='->',
    #    minlength=.1)
        
    plt.gca().set_aspect('equal')        
    plt.title('Velocity magnitude and forces')
    plt.xlabel('$r/r_g$')
    plt.ylabel('$r/r_g$')
    #plt.savefig(fileout)
    plt.show()

    
    print('Index:' + index_str)