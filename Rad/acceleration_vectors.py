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
                       

    
# Create meshgrid (in units of gravitational radii, then converted to physical)
ymin = 0
ymax = 100
xmin = 0
xmax = 100

grid_y, grid_x = np.mgrid[ymin:ymax:200j,xmin:xmax:200j]
grid_y = grid_y
grid_x = grid_x

## Loop through all simulation files
for index in range(1501,1502):
    if index < 1000:
        index_str = str('0' + str(index))
    else:
        index_str = str(index)
    
    #File operations
    filename = 'sim' + index_str
    #filein = open('/Users/Anton/Desktop/Data/d300a0/simavg0100-0189_rel.dat','rb')
    filein = open('/Users/Anton/Desktop/Data/d300a0/' + filename + '.dat','rb')
    #datafile_all_forces = '/Users/Anton/Desktop/Data/Binaries/RAD-600-d300a0/' + index_str + '.npy'
    #datafile_all_forces = '/Users/Anton/Desktop/Data/Binaries/RAD-avg-d300a0/average.npy'


    #fileout = '/Users/Anton/Desktop/Data/Image/' + filename + '.png'
    
    datain = np.loadtxt(filein,input_dtype.dtype_rel_radiation())
    filein.close()
    
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
    
     ## Velocities ##
          
    u_radial = (datain['u_1']/datain['u_t'])
    u_theta = (datain['u_2']/datain['u_t'])
    u_phi = (datain['u_3']/datain['u_t'])
    
    
   
    #grid_magnetic_pressure = griddata(points, datain['bsq']/2., (grid_x, grid_y), method='linear')
    #grid_gas_pressure = griddata(points, gas_pressure, (grid_x, grid_y), method='linear')
    grid_u_magnitude = griddata(points, np.sqrt(u_radial*u_radial + u_theta*u_theta*g22*g22), (grid_x, grid_y), method='linear')
    ux = (u_theta*np.cos(theta)*r + u_radial*np.sin(theta))
    uy = (-u_theta*np.sin(theta)*r + u_radial*np.cos(theta))
    grid_ux  = griddata(points, ux, (grid_x, grid_y), method='linear')
    grid_uy  = griddata(points, uy, (grid_x, grid_y), method='linear')
        
    fig = plt.figure()
    ax = fig.gca()
    plt.rc('text', usetex=True)
    
    ymin = 0
    ymax = 100
    xmin = 0
    xmax = 70
    
    plt.xlim(xmin=xmin, xmax=xmax)
    plt.ylim(ymin=ymin, ymax=ymax)
    
    c = plt.contourf(grid_x, grid_y, np.log10(grid_rho), extend='both', levels=np.linspace(-25,-19,61), cmap='gray', alpha=0.5)
    #cax, args = matplotlib.colorbar.make_axes(ax, location='left', pad=0.125)
    bar = fig.colorbar(c)#, cax=cax)
    #bar.ax.yaxis.set_ticks_position('left')
    #plt.sca(ax)
    
    lw = 1+grid_u_magnitude*10.
    s = plt.streamplot(grid_x, grid_y, grid_ux, grid_uy, 
        density=1, color='#444444', 
        arrowsize=2.5, arrowstyle='->',
        minlength=.1, linewidth=lw)
        
    ax.set_aspect('equal')
    #plt.title('Density and velocites')
    plt.xlabel('$r/r_g$')
    #plt.ylabel('$r/r_g$')
    plt.tick_params(axis='both', which='both', bottom='on', top='off', labelbottom='on', right='off', left='off', labelleft='off')

    fileout = '/Users/Anton/Dropbox/Aleksander/Figures/simavg0100-0189/density_velocities_'+filename+'.png'
    #fileout = '/Users/Anton/Dropbox/Aleksander/Figures/simavg0100-0189/density_velocities.png'
    plt.savefig(fileout, bbox_inches='tight')#, bbox_extra_artists=[cax]) 
    
    fig.show()

    
    print('Index:' + index_str)