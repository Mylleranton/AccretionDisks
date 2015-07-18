import pylab
import matplotlib.pyplot as plt
import matplotlib
import os
import numpy as np
from scipy.interpolate import griddata

from matplotlib.ticker import FormatStrFormatter

ifxlabel=0
ifylabel=0
ifcbar=0

##general preamble
matplotlib.rcdefaults()
matplotlib.rc('font',**{'family':'serif','size':16})
os.environ['PATH'] = os.environ['PATH'] + ':/opt/local/bin'
os.environ['PATH']

##files
#name='r001'
#number='0012-0018'

#filein1 = '/Users/asadowski/Sci/Koral/150614-3d-sims/'+name+'/simphisli0'+number+'.dat'
#filein1 = '/Users/asadowski/Sci/Koral/150614-3d-sims/r001/simavgphiavg0012-0018.dat'
#filein2 = '/Users/asadowski/Sci/Koral/150614-3d-sims/r003/simavgphiavg0012-0017.dat'
#filein3 = '/Users/asadowski/Sci/Koral/150614-3d-sims/r011/simavgphiavg0010-0015.dat'
#filein4 = '/Users/asadowski/Sci/Koral/150614-3d-sims/r020/simavgphiavg0012-0018.dat'
filein = '/Users/Anton/Desktop/Data/simavg0070-0134.dat'
fileout = '/Users/Anton/Desktop/Data/Analysis/image_001.png'

datain = pylab.loadtxt(filein)
        
##plotting preamble
plt.close('all')
plt.rc('text', usetex=True)
fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1,5, sharex=False, sharey=True, figsize=(14,6.5))
fig.subplots_adjust(bottom=0.1, left=0.05, right=0.92, top=0.9,wspace=0)

#########################
######## plot #1 ########
#########################


##plots
for i in range(1,6):
    
    if i==1:
        data=datain
        ax=ax1
        ax.text(15, 53.5,r'$\dot m=10$',fontsize=19, ha='center', va='top')
    if i==2:
        data=data2
        ax=ax2
        ax.text(15, 53.5,r'$\dot m=150$',fontsize=19, ha='center', va='top')
    if i==3:
        data=data3
        ax=ax3
        ax.text(15, 53,r'$\dot m=15\,\,a_*=0.7$',fontsize=19, ha='center', va='top')
    if i==4:
        data=data4
        ax=ax4
        ax.text(15, 53,r'$\dot m=10\,\,M_{\rm{BH}}=1000$',fontsize=19, ha='center', va='top')
    if i==5:
        data=data5
        ax=ax5
        ax.text(15, 53.5,r'$\dot m=10\,\,2D$',fontsize=19, ha='center', va='top')
        
    points=(data[:,4-1]*np.sin(data[:,5-1]),data[:,4-1]*np.cos(data[:,5-1]))
    radflux=np.log10(np.sqrt(data[:,15-1]*data[:,15-1] + data[:,16-1]*data[:,16-1]*data[:,4-1]*data[:,4-1]))
    
    if i==4:
        radflux=radflux+2.
    
    r=data[:,4-1]
    th=data[:,5-1]
    Fr=data[:,15-1]
    Fth=data[:,16-1]*r
    Fx=Fr*np.sin(th)+Fth*np.cos(th)
    Fy=Fr*np.cos(th)-Fth*np.sin(th)
     
    vminmy = 22.6
    vmaxmy = 26
    
    #delete the BH
    for j in range(0,len(r)):
        if r[j]<2.:        
            radflux[j]=vminmy
     
    xlim=30
    ylim1=-12.5
    ylim2=50
    
    
    grid_y, grid_x = np.mgrid[ylim1:ylim2:300j, 0:xlim:150j]
    
    grid_z2 = griddata(points, radflux, (grid_x, grid_y), method='linear')
    
    grid_Fx = griddata(points, Fx, (grid_x, grid_y), method='linear')
    grid_Fy = griddata(points, Fy, (grid_x, grid_y), method='linear')
    grid_Fr = griddata(points, Fr, (grid_x, grid_y), method='linear')
    
    
    
    
    im = ax.pcolor(grid_x, grid_y, grid_z2, 
        vmin=vminmy,vmax=vmaxmy,cmap=plt.cm.YlOrBr_r)
    ax.set_ylim((ylim1,ylim2))
    ax.set_xlim((0,xlim))
        
    lw = 4*(grid_z2-22)/4+0.5
    sim1 = ax.streamplot(grid_x, grid_y, -grid_Fx, -grid_Fy, density=.5, color='#dbe8e8', 
        arrowsize=2, arrowstyle='->',
        minlength=0.1,linewidth=lw)
        
    ax.set_xticks([0,10,20])
    
    if i==1:
        ax.set_ylabel(r'$z/M$',fontsize=22,labelpad=-10)
    if i==3:
        ax.set_xlabel(r'$r/M$',fontsize=22,labelpad=0)

    cntr=ax.contour(grid_x, grid_y, grid_Fr, [0],linewidths=4, colors = "#aba8e8",linestyles='--')
    

##axes

cbaxes = fig.add_axes([0.93, 0.1, 0.01, 0.8]) 
cbar = fig.colorbar(im,cax=cbaxes,format='%3.3g')
cbar.set_label(r'$\log\,F_{\rm rad}\,\,\rm[erg/s/cm^2]$',fontsize=20,labelpad=3)    
cbar.ax.tick_params(labelsize=12) 
#cbar.ax.yaxis.set_ticks_position('left')

##save
plt.savefig(fileout)