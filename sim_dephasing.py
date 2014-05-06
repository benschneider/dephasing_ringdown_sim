# -*- coding: utf-8 -*-
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

#import equations as eq
import parsers as parser
import simulation_functions as sim
'''
Note if debugging these files, better use execfile('...py') 
otherwise changes are updated only after a kernel restart
still missing some pieces of code
'''

# -------- Set variabes preperation for simulation ---------

#output filenames
filename1 = 'QsQr2.txt'
filename2 = 'qsfit2.mtx'
filename3 = 'qrfit2.mtx'
filename4 = 'sim_dephasing2.mtx'
 
#input filename
filename_lockinfilter = 'experimental/1_2_filt8.Wfm.csv' #100us long 1000pt res

#set variables
Qs = 1500
w_res = 300.0*2*np.pi

#Ringdowntime in usec
t_0 = -40 #start
t_1 = 40 #int(Qs/w_res*30.0) #stop
t_p = 801 #number of points

#Frequency in MHz
span = w_res/Qs*10
w_0 = w_res-span
w_1 = w_res+span
w_p = 1501

#Ringdown Q factor range
Qr_0 = 1500
Qr_1 = 10000
Qr_p = 17

#make np arrays
t_array  = np.linspace(t_0, t_1, t_p)
w_array  = np.linspace(w_0, w_1, w_p)
Qr_array = np.linspace(Qr_0, Qr_1, Qr_p)

# -------- Simulation and fitting ---------

#calculate required dephasing:
Qd_array = sim.eq.Qd(Qs,Qr_array[1:])

#create non dephased matrix:
matrix3d = sim.get_ringdown_matrix_Y(w_array, w_res, t_array, Qr_array) 

t1 = time()
print 'dephase matrix'
matrix3d = sim.get_dephased_matrix(matrix3d, Qd_array, w_array, method='lor') 
#setting method to gaus uses a different method and is a lot faster
print time()-t1

#grab lockin response and extract lockin filter function:
lockin_response = parser.loadcsv(filename_lockinfilter)
lockin_filter = sim.filter2d_lockin(lockin_response) #calculate and return a normalized filter
#the following is done by manual selection
#pl.plot(lockin_filter[207:207+t_p])

filt_fun = lockin_filter[0:450] #[207:207+t_p]
#filtfun = lockin_filter
#tmp = np.zeros([1,len(filtfun)]) #store filter in the right format
#for j in range(0, tmp.shape[0]):
#    tmp[j] = filtfun #store normalized filte
#
#filtfun = tmp

print 'convolve matrix with filter function'
matrix3d = sim.get_matrix_lockin_convolved(matrix3d, filt_fun)
#A = sim.convolution_1d(matrix3d[0],filt_fun)

'''
print 'fit ringdown/spectral Q for each point of dephasing'
Qrs, qrfit = sim.fit_mat_ringdown(t_array, matrix3d) 
Qsp, qsfit = sim.fit_matrix_spectral(w_array, matrix3d) 
'''

''' 
qrfit contains a 2d data with alternating lines 
containing the data and the fit 
Qrs contains the Q factor and dephasng 
'''


print ' ------- Save data into files --------'

'''
#save fiting results into dat file
stuff = (Qr_array, Qrs[0], Qsp[0])
parser.savedat(filename1, stuff, delimiter='\t') 

#Save the Spectral fitting matix into one MTX file
#header file
head = ['Units', 'Qs_fits',
        'RF frequency (MHz)', str(w_array[0]/2/np.pi), str(w_array[-1]/2/np.pi),
        'Ringd_Q', str(Qr_array[-1]) , str(Qr_array[0]),
        'none', '0', '1']
parser.savemtx(filename2, qsfit, header=head) #save in MTX format

#Save the Ringdown fitting matix into one MTX file
head = ['Units', 'Qs_fits',
        'Time (us)', str(t_array[0]), str(t_array[-1]),
        'Ringd_Q', str(Qr_array[-1]) , str(Qr_array[0]),
        'none', '0', '1']
parser.savemtx(filename3, qrfit, header=head)
'''

#Save the dephasing matix as MTX file
head = ['Units', 'X Y R [V]#values',
        'Time (us)', str(t_array[0]), str(t_array[-1]),
        'RF frequency (MHz)', str(w_array[-1]/2/np.pi), str(w_array[0]/2/np.pi),
        'Ringd_Q', str(Qr_array[0]), str(Qr_array[-1])]
parser.savemtx(filename4, matrix3d, header = head)



'''
#extract the lockin convolution filter from (a step function)
mtx1_data, mtx1_header = mp.loadmtx('meas_data/n20dBm_amp.mtx')
mtx1_data = mtx1_data[0] #it is a 2d colour map no 3d req.
data_off = mtx1_data[0] - mtx1_data[0].min()
data_off = data_off/data_off.max()
data_filt = -np.diff(data_off,1) #obtain Filter
data_filt2 = data_filt[2000:4000] #croped for this particular data set

#averaging and croping some parts of a filter
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