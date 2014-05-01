'''
Ben Schneider

this script generates a simulated ringdown dataset
be carefull not to increase the limits too much a high
(in consumes space and memory)
'''

#from numpy import sin, cos, e, pi, sqrt, linspace, array, zeros, ones, savetxt
import numpy as np
#from time import time
#from scipy.optimize import curve_fit
#from scipy import signal
#import scipy.signal as signal
#from scipy.signal.convolve
#import scipy.ndimage #The fastest filter
import matplotlib.pyplot as pl #for debugging

import equations as eq
import mtxparser as mtx
import simulation_functions as sim
'''Note if debugging these files, better use execfile('...py') 
otherwise changes are updated only after a kernel restart
'''

filename = 'sim_dephasing.mtx'

#set variables
Q = 6000.0 #a first assumption of the Q
Qs = 1500
w0 = 300.0*2*np.pi
dephasing_points = 3

dephasing = np.linspace(0, 1, dephasing_points)#in (Mhz*2.355) (FWHM)
span = w0/Qs*10 
rdtime = int(Q/w0*12.0)

t = np.linspace(-3, rdtime, 5001) #Ringdowntime in usec
w = np.linspace(w0-span, w0+span, 1501) #Frequency range
#matrix = np.zeros((dephasing_points, len(w), len(t))) #The result matrix


'''
#set some range of Qr between:
#Qs = 1500
#Qr = Qs and Qd = 0
#to Qr = 10*Qs and Qs = ...
#start with:
'''

Qr = np.linspace(Qs,Qs*10,dephasing_points)
dephasing = eq.Qd(Qs,Qr)

matrix = sim.get_ringdown_Y2(w, w0, t, Qr)

#matrix[0] = sim.get_ringdown_Y(w, w0, t, Q) #get non dephased matrix

#matrix = sim.get_dephased_matrix(matrix, dephasing) #get dephasing to matrix
matrix = sim.get_dephased_matrix_Y(matrix, dephasing, w)
#fit the spectral and ringdown Q of the matrix 
'''
qrfit contains a 2d data with alternating lines containing the data and the fit
Qrs contains the Q factor and dephasng
'''
Qrs, qrfit = sim.fit_mat_ringdown(t, matrix) #fit ringdown Q for each dephasing
Qsp, qsfit = sim.fit_matrix_spectral(w, matrix) #fir spectral Q for each dephasing



#Dephasing (per pixel) per (w0/Q)
#deph_fact = ((w.max()-w.min())/len(w))*Q/w0
deph_fact = 1
stuff = (dephasing*deph_fact, Qrs[0], Qsp[0])
mtx.savedat("QsQr.txt", stuff, delimiter='\t') #save fiting results into dat file

#Save the Spectral fitting matix into one MTX file
head = ['Units', 'Qs_fits',
        'Dephasing (rel to FWHM)', '0', str(deph_fact*dephasing[-1]),
        'RF frequency (MHz)', str(w[-1]/2/np.pi), str(w[0]/2/np.pi),
        'other', '0', '1']

mtx.savemtx('qsfit.mtx', qsfit, header=head) #save as MTX

#Save the Ringdown fitting matix into one MTX file
head = ['Units', 'Qs_fits',
        'Dephasing (rel to FWHM)', '0', str(deph_fact*dephasing[-1]),
        'Time (us)', str(t[-1]), str(t[0]),
        'other', '0', '1']
mtx.savemtx('qrfit.mtx', qrfit, header=head) #save as MTX

#Save the dephasing matix into another MTX file
head = ['Units', 'X Y R [V]#values',
        'Time (us)', str(t[0]), str(t[-1]),
        'RF frequency (MHz)', str(w[-1]/2/np.pi), str(w[0]/2/np.pi),
        'Dephasing (rel to FWHM)', '0', str(deph_fact*dephasing[-1])]
mtx.savemtx(filename, matrix, header = head) #save as MTX


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