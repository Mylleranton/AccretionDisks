import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mline
import matplotlib
import os

matplotlib.rcdefaults()
matplotlib.rc('font',**{'family':'serif','size':16})
os.environ['PATH'] = os.environ['PATH'] + ':/opt/local/bin'
os.environ['PATH']


datafile = '/Users/Anton/Desktop/Data/Binaries/RAD-avg-d300a0/average.npy'
savefile = '/Users/Anton/Desktop/Data/Binaries/RAD-avg-d300a0/average_plot.npy'

fileout = '/Users/Anton/Desktop/Data/Image/average-RAD-d300a0_3.png'
index = 1
try:
    loaded_data = np.load(datafile)
except IOError as error:
    print('Error in loading datafile:', error)
    
loaded_array = loaded_data[:,2:14]
coordinates = loaded_data[:,0:2]



#for i in range(1002,1013):
#    if i < 1000:
#        index_str = str('0' + str(i))
#    else:
#        index_str = str(i)
#        
#    index = index + 1
#    datafile = '/Users/Anton/Desktop/Data/Binaries/NO-RAD-600-d300a0/' + index_str + '.npy'
#    try:
#        loaded_arr = np.load(datafile)
#    except IOError as error:
#        print('Error in loading datafile:', error)
#    
#    loaded_array = np.add(loaded_array, loaded_arr[:,2:14])

loaded_array = np.divide(loaded_array, index)
loaded_array2 = np.concatenate((coordinates,loaded_array),axis=1)
np.save(savefile, loaded_array2)

#########
######### Plot
#########

# Wind = (40,70), Jet = (15,80), Disk = (80,20)
x         = np.array([15,45,60])
X_all     = np.array([x[0],x[0],x[0],x[0],x[0],x[0], x[1],x[1],x[1],x[1],x[1],x[1], x[2],x[2],x[2],x[2],x[2],x[2]])
y         = np.array([40,55,20])
Y_all     = np.array([y[0],y[0],y[0],y[0],y[0],y[0], y[1],y[1],y[1],y[1],y[1],y[1], y[2],y[2],y[2],y[2],y[2],y[2]])

row_1 = loaded_array2[np.logical_and(loaded_array2[:,0]==x[0], loaded_array2[:,1]==y[0])]
row_2 = loaded_array2[np.logical_and(loaded_array2[:,0]==x[1], loaded_array2[:,1]==y[1])]
row_3 = loaded_array2[np.logical_and(loaded_array2[:,0]==x[2], loaded_array2[:,1]==y[2])]

gradients_all=np.array([[row_1[0,2], row_1[0,3]],
                        [row_1[0,4], row_1[0,5]],
                        [row_1[0,6], row_1[0,7]],
                        [row_1[0,8], row_1[0,9]],
                        [row_1[0,10], row_1[0,11]],
                        [row_1[0,12], row_1[0,13]],
                        [row_2[0,2], row_2[0,3]],
                        [row_2[0,4], row_2[0,5]],
                        [row_2[0,6], row_2[0,7]],
                        [row_2[0,8], row_2[0,9]],
                        [row_2[0,10], row_2[0,11]],
                        [row_2[0,12], row_2[0,13]],
                        [row_3[0,2], row_3[0,3]],
                        [row_3[0,4], row_3[0,5]],
                        [row_3[0,6], row_3[0,7]],
                        [row_3[0,8], row_3[0,9]],
                        [row_3[0,10], row_3[0,11]],
                        [row_3[0,12], row_3[0,13]]])



## Plot the interpolated data

plt.figure()
plt.rc('text', usetex=True)

ymin = 0
ymax = 100
xmin = 0
xmax = 100

plt.xlim(xmin=xmin, xmax=xmax)
plt.ylim(ymin=ymin, ymax=ymax)


q_all = plt.quiver(X_all, Y_all, gradients_all[:,0], gradients_all[:,1], 
    color=['blue', 'green', 'yellow', 'red','cyan', 'black'], 
    width=.004, scale=25)


p1 = mline.Line2D([], [], color='blue', label='Gravity')
p2 = mline.Line2D([], [], color='green', label='Gas')
p3 = mline.Line2D([], [], color='yellow', label='Magnetic pressure')
p4 = mline.Line2D([], [], color='red', label='Centrifugal')
p5 = mline.Line2D([], [], color='cyan', label='Radiation')
p6 = mline.Line2D([], [], color='black', label='Total')

plt.legend(handles=[p1,p2,p3,p4,p5,p6], fontsize='x-small')

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