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
for index in range(1120,1121):
    if index < 1000:
        index_str = str('0' + str(index))
    else:
        index_str = str(index)
    
    #File operations
    filename = 'simext' + index_str
    #filein = open('/Volumes/Seagate/4Anton/d300a0/dt100/simavg0100-0189_rel.dat','rb')
    filein = open('/Volumes/Seagate/4Anton/d300a0/dt10/'+filename+'.dat','rb')
    #datafile_all_forces = '/Users/Anton/Desktop/Data/Binaries/RAD-600-d300a0/' + index_str + '.npy'
    #datafile_all_forces = '/Users/Anton/Desktop/Data/Binaries/RAD-avg-d300a0/average.npy'


    #fileout = '/Users/Anton/Desktop/Data/Image/' + filename + '.png'
    
    datain = np.loadtxt(filein,input_dtype.dtype_rel_radiation())
    filein.close()
    
    #Convert coordinates from spherical to cartesian (physical units)
    r = datain['r']
    theta = datain['theta']
    phi = datain['phi']
    ehat = datain['rad_ehat']
    
    points = (r*np.sin(theta),r*np.cos(theta))
    grid_ehat = griddata(points, ehat, (grid_x, grid_y), method='linear')
    
    ## Metric determinant ##
    a=0
    g_sch = np.power(r, 2)*np.sin(theta)
    
    g00 = -(1-np.divide(2.*r,(np.power(r, 2)+a**(2.)*np.power(np.cos(theta),2))))
    g11 = np.sqrt((np.divide( (np.power(r, 2)+a**(2.)*np.power(np.cos(theta),2)), (np.power(r, 2)-2.*r+a**2) ) ))
    g22 = np.sqrt((np.power(r,2)+a**(2.)*np.power(np.cos(theta),2)))
    g33 = np.sqrt(np.power(np.sin(theta), 2)*(   np.power(r, 2)  + a**(2.) + np.divide(2.*r, (np.power(r, 2)-2.*r+a**(2.)))*a**(2.)*np.power(np.sin(theta),2)))
    
    ## Radflux ##
    rad_radial = -datain['R_10']*g11
    rad_theta = -datain['R_20']*g22
    
    rad_x = (rad_theta*np.cos(theta) + rad_radial*np.sin(theta))
    rad_y = (-rad_theta*np.sin(theta) + rad_radial*np.cos(theta))
    grid_rad_x  = griddata(points, rad_x, (grid_x, grid_y), method='linear')
    grid_rad_y  = griddata(points, rad_y, (grid_x, grid_y), method='linear')
        
    grid_rad_magnitude = griddata(points, np.sqrt(rad_radial*rad_radial + rad_theta*rad_theta), (grid_x, grid_y), method='linear')

    
    fig = plt.figure()
    ax = fig.gca()
    plt.rc('text', usetex=True)
    
    ymin = 0
    ymax = 100
    xmin = 0
    xmax = 70
    
    plt.xlim(xmin=xmin, xmax=xmax)
    plt.ylim(ymin=ymin, ymax=ymax)
    
    c = plt.contourf(grid_x, grid_y, np.log10(grid_rad_magnitude*1.66e47), extend='both', levels=np.linspace(21,25,61), cmap='gist_heat', alpha=0.5)
    cax, args = matplotlib.colorbar.make_axes(ax, location='top', shrink=0.45)
    bar = fig.colorbar(c, cax=cax, ticks=np.linspace(21,25,5), orientation='horizontal')
    bar.ax.xaxis.set_ticks_position('top')    
    plt.sca(ax)
    
    lw = np.abs(np.log10(grid_rad_magnitude))-20

    s = plt.streamplot(grid_x, grid_y, grid_rad_x, grid_rad_y, 
        density=1, color='#222222', 
        arrowsize=2.5, arrowstyle='->',
        minlength=.1, linewidth=2)
        
    ax.set_aspect('equal')
    #plt.title('Density and velocites')
    plt.xlabel('$x/r_g$')
    #plt.ylabel('$r/r_g$')
    plt.tick_params(axis='both', which='both', bottom='on', top='off', labelbottom='on', right='off', left='off', labelleft='off')

    fileout = '/Users/Anton/Dropbox/Aleksander/Figures/simavg0100-0189/radiation_flux_'+filename+'.png'
    #fileout = '/Users/Anton/Dropbox/Aleksander/Figures/simavg0100-0189/density_velocities.png'
    fig.savefig(fileout, bbox_inches='tight') 
    
    fig.show()

    
    print('Index:' + index_str)