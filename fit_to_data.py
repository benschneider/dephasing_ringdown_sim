# -*- coding: utf-8 -*-
"""
Created on Tue May  6 17:54:15 2014

@author: benschneider

structure of this code:
1st.
load simulated data with the fixed Qs and different Qrs. 
(simulated data set needs to have the same time/point ratio as the measured data)
load measured data with known Qs and unknown Qr 

2nd.
for each (Qr) in simulated data and the measured data list:
normalize (amplitude from 0 to 1) for measured and simulated data.

3rd.
find position in time where amplutude is 0.5
define this position to be at t = 0

4th.
compare (or fit), simulated with measured line for each Qr

5th.
show/print result
"""

import parsers as parser
import matplotlib.pyplot as pl
import numpy as np

def norm_line(linecut):
    '''set the number range from 0 to 1'''
    tmp = linecut - linecut.min()
    linecut = tmp/tmp.max()
    return linecut

def find_index(linecut,target):
    '''go through each element and 
    see if it is closer to the target element.
    if so, store position, 
    otherwise continue searching until all elements have been checked.    
    '''    
    i = 0
    tmp1 = 1
    for num in linecut:
        value = abs(target-num)
        if value < tmp1:
            tmp1 = value
            index = i
        i+=1
    
    return index#,tmp1


filename1 = 'sim_dephasing2.mtx' #simulated data (3d np array)
filename2 = 'experimental/n20dBm_amp_subed.l.139.linecut.dat'
#filename2 = 'experimental/n20dBm_amp_for_fitting.l.139.linecut.dat'

meas_raw = parser.loaddat(filename2) #load data
sim_raw,  sim_raw_head  = parser.loadmtx(filename1) #load data

fit_adj = 0.5
fit_adjf = 0.99
fit_range = 110 #points
fit_left = 200 #points
fit_right = 250 #points

#----- Measured trace ---- (no loop required)
meas_data = np.array(meas_raw[1]) 
meas_data = norm_line(meas_data) #normalize
meas_time = np.array(meas_raw[0])#time in sec
meas_time = meas_time*1e6 #set time in usec

#find zero position and adjust time
meas_data = norm_line(meas_data) #pre normalize data
meas_0pos1 = find_index(meas_data, fit_adj)
meas_0pos = find_index(meas_data[meas_0pos1-60:meas_0pos1+100], fit_adjf)
meas_0pos = meas_0pos+meas_0pos1-60
meas_0off = meas_time[meas_0pos]
meas_time = (meas_time - meas_0off)

#crop around zero postion
meas_p0 = (meas_0pos - fit_left)
meas_p1 = (meas_0pos + fit_right)
meas_time = meas_time[meas_p0:meas_p1]
meas_data = meas_data[meas_p0:meas_p1]

#post normalize data
meas_data = norm_line(meas_data)

pl.figure(1)
pl.close()
pl.figure(2)
pl.close()

pl.figure(1)
pl.plot(-meas_data+1)


#----- Simulated trace ---- (loop is required)

sim_res_index = round(sim_raw.shape[1]/2) #3d mat resonance pos

#produce time axis for simulated data set
sim_t_0 = eval(sim_raw_head[3]) #start in usec
sim_t_1 = eval(sim_raw_head[4])
sim_time2 = np.linspace(sim_t_0, sim_t_1, sim_raw.shape[2])

Ki2 = np.zeros(sim_raw.shape[0])

for Qr_index in range(0,sim_raw.shape[0]):
    #Qr_index = 0 #this is temporarily
    sim_data = sim_raw[Qr_index][sim_res_index] #select Qr trace
    sim_data = sim_data[120:-1] #precrop for bugfixing
    sim_time = sim_time2[120:-1]
    
    #find zero position and adjust time
    sim_data = norm_line(sim_data) #pre normalize data
    sim_0pos1 = find_index(sim_data, fit_adj)
    sim_0pos = find_index(sim_data[sim_0pos1-60:sim_0pos1+100], fit_adjf)
    sim_0pos = sim_0pos+sim_0pos1-60

    sim_0off = sim_time[sim_0pos]
    sim_time = (sim_time - sim_0off)
    
    #crop around zero postion
    sim_p0 = (sim_0pos - fit_left) 
    sim_p1 = (sim_0pos + fit_right)
    sim_data = sim_data[sim_p0:sim_p1]
    sim_time = sim_time[sim_p0:sim_p1]
    
    #post normalize data
    #sim_data = norm_line(sim_data)
    
    #calculate difference
    Ki2_tmp= (meas_data-sim_data)**2
    Ki2[Qr_index] = Ki2_tmp.sum()

    #pl.close()
    #pl.plot(sim_time,sim_data)
    #pl.plot(meas_time,meas_data)
    pl.figure(1)
    pl.plot(-sim_data+1)


pl.figure(2)
pl.plot(Ki2)