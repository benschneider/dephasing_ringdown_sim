'''
Ben Schneider

this script generates a simulated ringdown dataset
be carefull not to increase the limits too much a high
(in consumes space and memory)
'''

#from numpy import sin, cos, e, pi, sqrt, linspace, array, zeros, ones, savetxt
import numpy as np
#from scipy.optimize import curve_fit
#from scipy import signal
#import scipy.signal as signal
#from scipy.signal.convolve
#import scipy.ndimage #The fastest filter

#for debugging
from time import time
import matplotlib.pyplot as pl 

import equations as eq
import mtxparser as mtx
import simulation_functions as sim
'''Note if debugging these files, better use execfile('...py') 
otherwise changes are updated only after a kernel restart
still missing some pieces of code
'''

# -------- Set variabes preperation for simulation ---------

#output filenames
filename1 = 'QsQr.txt'
filename2 = 'qsfit.mtx'
filename3 = 'qsfit.mtx'
filename4 = 'sim_dephasing.mtx'
 
#set variables
Qs = 1500
w0 = 300.0*2*np.pi

#Ringdowntime in usec
t_0 = -3 #start
t_1 = int(Qs/w0*20.0) #stop
t_p = 3001 #number of points

#Frequency in MHz
span = w0/Qs*10
w_0 = w0-span
w_1 = w0+span
w_p = 1501

#Ringdown Q factor range
Qr_0 = Qs
Qr_1 = Qs*30
Qr_p = 3

#make np arrays
t  = np.linspace(t_0, t_1, t_p)
w  = np.linspace(w_0, w_1, w_p)
Qr = np.linspace(Qr_0, Qr_1, Qr_p)

# -------- Start simulation and fitting ---------

#calculate required dephasing
Qd = eq.Qd(Qs,Qr[1:]) #get array

#create non dephased matrix
matrix = sim.get_ringdown_Y2(w, w0, t, Qr) 

#normalize matrix
matrix = sim.get_normed_matrix(matrix, Qr)

#dephase matrix accordingly
t1 = time()
matrix = sim.get_dephased_matrix(matrix, Qd, w, method='lor')
print time()-t1

#fit ringdown Q for each point of dephasing 
Qrs, qrfit = sim.fit_mat_ringdown(t, matrix) 
#fit spectral Q for each point of dephasing
Qsp, qsfit = sim.fit_matrix_spectral(w, matrix) 

''' 
qrfit contains a 2d data with alternating lines 
containing the data and the fit 
Qrs contains the Q factor and dephasng 
'''


# ------- Save data into files --------


#save fiting results into dat file
stuff = (Qr, Qrs[0], Qsp[0])
mtx.savedat(filename1, stuff, delimiter='\t') 

#Save the Spectral fitting matix into one MTX file
#header file
head = ['Units', 'Qs_fits',
        'Ringd_Q', str(Qr[0]) , str(Qr[-1]),
        'RF frequency (MHz)', str(w[-1]/2/np.pi), str(w[0]/2/np.pi),
        'none', '0', '1']
mtx.savemtx(filename2, qsfit, header=head) #save in MTX format

#Save the Ringdown fitting matix into one MTX file
head = ['Units', 'Qs_fits',
        'Ringd_Q', str(Qr[0]) , str(Qr[-1]),
        'Time (us)', str(t[-1]), str(t[0]),
        'none', '0', '1']
mtx.savemtx(filename3, qrfit, header=head)

#Save the dephasing matix into another MTX file
head = ['Units', 'X Y R [V]#values',
        'Time (us)', str(t[0]), str(t[-1]),
        'RF frequency (MHz)', str(w[-1]/2/np.pi), str(w[0]/2/np.pi),
        'Ringd_Q', str(Qr[0]), str(Qr[-1])]
mtx.savemtx(filename4, matrix, header = head)



'''
#extract the lockin convolution filter from (a step function)
mtx1_data, mtx1_header = mp.loadmtx('meas_data/n20dBm_amp.mtx')
mtx1_data = mtx1_data[0] #it is a 2d colour map no 3d req.
data_off = mtx1_data[0] - mtx1_data[0].min()
data_off = data_off/data_off.max()
data_filt = -np.diff(data_off,1) #obtain Filter
data_filt2 = data_filt[2000:4000] #croped for this particular data set

#averaging and croping some parts of the filter
fact = 1
a = np.zeros(len(t)/fact)
a = np.array(a)
avg = float(len(data_filt))/float(len(t))*fact
avg = int(avg)
j = 0
j2 = 0
for i in range(0,len(data_filt)+1):
    j = i%avg
    if j == 0 and i != 0:
        j2 +=1
    if (j+j2) < len(t)/fact:
        a[j+j2] += data_filt[i]/avg
    else:
        break;
data_filt3 = a[170:330] #crop an averaged filter
'''