# -*- coding: utf-8 -*-
'''
Ben Schneider

this script generates a simulated ringdown dataset
be carefull not to increase the limits too much.
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
#from time import time
#import matplotlib.pyplot as pl 

#import equations as eq
import parsers as parser
import simulation_functions as sim

execfile('parameters.py') #this file contain all the parameters above

#make np arrays
t_array  = np.linspace(t_0, t_1, t_p)
w_array  = np.linspace(w_0, w_1, w_p)
Qr_array = np.linspace(Qr_0, Qr_1, Qr_p)
Qd_array = sim.eq.Qd(Qs,Qr_array[1:]) #dephasing matrix

print 'calculate base matrix'
matrix3d = sim.get_ringdown_matrix_Y(w_array, w_res, t_array, Qr_array) 

print 'dephase matrix'
matrix3d = sim.get_dephased_matrix(matrix3d, Qd_array, w_array, w_res, method='lor') 

print 'extract lockin filter' 
#grab lockin response and extract lockin filter function:
'''
lockin_response = parser.loadcsv(filename_lockinfilter)
lockin_response = lockin_response[0]
'''
lockin_response = parser.loaddat(filename_3)
lockin_val = np.array(lockin_response[1])

#calculate and return a normalized and shrinked filter
lockin_filter = sim.filter2d_lockin(lockin_val) 
#lockin_filter = sim.shrink_extend_array(lockin_filter, t_p)

#the following cropping of the filter, is done by manual selection
filt_peakidx = lockin_filter.argmax()
filt_crop = (filt_peakidx*2)
filt_fun = lockin_filter[0:filt_crop]

#renormalize croped filter
filt_fun = sim.norm_filter(filt_fun) 

print 'convolve matrix with filter function'
matrix3d = sim.get_matrix_lockin_convolved(matrix3d, filt_fun)


print ' ------- Save data as MTX file --------'
head = ['Units', 'X Y R [V]#values',
        'Time (us)', str(t_array[0]), str(t_array[-1]),
        'RF frequency (MHz)', str(w_array[0]/2/np.pi), str(w_array[-1]/2/np.pi), 
        'Ringd_Q', str(Qr_array[0]), str(Qr_array[-1])]
parser.savemtx(filename_1, matrix3d, header = head)