# -*- coding: utf-8 -*-
"""
Created on Tue Apr  8 21:04:32 2014

@author: benschneider
"""
#from scipy import signal
import numpy as np
#from time import time
#from scipy.optimize import curve_fit
import scipy.signal as signal
#import scipy.ndimage as ndimage #A fast filter
#import pylab as pl
import matplotlib
#matplotlib.use('pdf')
import matplotlib.pylab as pl
#from os import system as terminal
%matplotlib inline

execfile('mtx_parser.py')
execfile('equations.py') #load equations in

''' #for loading two individual files
import laddat #to simply load or save dat files
file1 = 'n20dBm_amp.l.0.linecut.dat'
file2 = 'n20dBm_amp.l.140.linecut.dat'
class resonance_data(object):
    off_res = np.array(laddat(file1)[1])
    off_time = np.array(laddat(file1)[0])
    on_res = np.array(laddat(file2)[1])
    on_time = np.array(laddat(file2)[0])
'''

mtx1_data, mtx1_header = loadmtx('n20dBm_amp.mtx')
mtx1_data = mtx1_data[0] #it is a 2d colour map so no 3d array required

#show fig
def show_f(filename1 = 'temp.pdf'):
    pl.savefig(filename1) #save as pdf
    #terminal('open '+filename1)

#from spyview I know line 0 is off resonant and line 130 is on resonant
#plotting those here confirms this
pl.figure(1)
pl.hold(True)
pl.subplot(2, 1, 1);
pl.plot(mtx1_data[0])
pl.plot(mtx1_data[140]) #this has a lower amplitude
pl.plot((mtx1_data[140]-mtx1_data[0]))
pl.subplot(2, 1, 2);
pl.imshow(mtx1_data)
pl.hold(False)
show_f('ON-OF-Diff.pdf')
pl.close()

data_on = mtx1_data[140] - mtx1_data[140].min()
data_on = data_on / data_on.max() -0.03
data_on = data_on / data_on.max()
data_off = mtx1_data[0] - mtx1_data[0].min()
data_off = data_off/data_off.max()
data_diff = mtx1_data[140] - mtx1_data[0]
data_diff = data_diff - data_diff.min()
data_diff = data_diff / data_diff.max()
data_t = np.linspace(eval(mtx1_header[3]),eval(mtx1_header[4]),len(mtx1_data[0]))
data_1 = np.hstack((data_on[::-1],data_on))
data_2 = np.hstack((data_off[::-1],data_off))
data_filt = -np.diff(data_off,1) #obtain Filter
data_filt2 = data_filt[2500:] #compressed version
data_range = int(len(data_off)/2)#-2500

#make a step function
def step(t,pos):
    j = np.zeros(len(t))
    for i in t:        
        if i < pos:
            j[i] = float(1.0)
        else:
            j[i] = float(0)
    return j

def mech(t,pos,tau):
    j = np.zeros(len(t))
    for i in t:        
        if i < pos:
            j[i] = float(1.0)
        else:
            j[i] = np.exp(-(float(i)-pos)/(2*tau))
    return j

st_t = np.array(range(0,len(data_off)))
st_p = len(st_t)/2

#make step function
data_step = step(st_t,st_p)
data_step = np.array(data_step) #stepfunction is ready
data_step2 = np.hstack((-data_step[::-1],data_step))

#make a mech decay simulation
tau = 111#1e-99
data_mech = mech(st_t,st_p,tau)
data_mech = np.array(data_mech)

#plot mech decay simulation and its filtured version
data_conv2 = signal.convolve(data_mech,data_filt)
data_conv2 = np.array(data_conv2[data_range:-data_range+2])
data_conv2 = data_conv2 - data_conv2.min()+0.03
data_conv2 = data_conv2 / data_conv2.max()

pl.figure(2)
#pl.plot(data_t*1e6,-data_diff+1, label='On res')
#pl.plot(data_t*1e6,data_conv2, label='conv')
pl.hold(True)
pl.plot(data_t[1500:4500]*1e6,-data_diff[1500:4500]+1, label='On res')
pl.plot(data_t[1500:4500]*1e6,data_conv2[1500:4500], label='conv')
#pl.plot(data_t*1e6,data_mech)
show_f('ringdown.pdf') #save and open
pl.hold(False)
pl.close()

pl.figure(3)
pl.hold(True)
pl.plot(data_filt)
pl.hold(False)
show_f('filter.pdf')
pl.close()

'''
#do a simple test of the filter:
data_conv = signal.convolve(data_step,data_filt2)
data_conv = np.array(data_conv)
data_conv = data_conv - data_conv.min()
data_conv = data_conv/data_conv.max()

#plot step convolution with filter and measured filtered data
#pl.figure(1)
#pl.plot(data_conv)
#pl.plot(data_off[2000:]/data_off[2000:].max())
'''



#test the filter
#data_deconv = signal.deconvolve(data_off,data_filt2)
#remove nan
'''
data_deconv = signal.deconvolve(data_3, data_filt2)
data_deconv2 = data_deconv[0]
for i in range(0,len(data_deconv[0])):
    if isnan(data_deconv[0][i]) == True:
        data_deconv2[i] = 0.0        
    elif isinf(data_deconv[0][i]) == True:
        data_deconv2[i] = 1.0e308
data_deconv2 = np.array(data_deconv2)
data_deconv2 = data_deconv2/data_deconv2.max()
'''