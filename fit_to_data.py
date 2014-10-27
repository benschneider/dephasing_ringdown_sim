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
#import simulation_functions as sim

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

def find_idx2(data, target, pre_target = 0.5, pre_range = 60):
    '''specifically for this type of ringdown data'''
    #pre normalize data
    data = norm_line(data) 
    
    #find pre position    
    pre_pos0 = find_index(data, pre_target)
    
    #select target area
    pre_0 = (pre_pos0 - pre_range)
    if pre_0 < 0: print 'pre_0 is less than 0, decrease pre_range value'
    pre_1 = (pre_pos0 + pre_range)
    data1 = data[pre_0:pre_1]

    #find target in target area
    pre_pos1 = find_index(data1, target) 
    pos = pre_0 +pre_pos1
    
    return pos
    
def crop_at_target(data1d, pos, fit_left, fit_right):
    p0 = pos - fit_left
    p1 = pos + fit_right
    if p0 < 0: 
        print 'p0 is less than 0, decrease left range'
    data1d2 = data1d[p0:p1]    
    return data1d2


#all manual changes are done in parameters.py
execfile('parameters.py') #this file contain all the parameters above


#----- Load files (given by parameters.py)----
meas_raw = parser.loaddat(filename_2)
sim_raw,  sim_raw_head  = parser.loadmtx(filename_1)


#to be adjusted to ensure the pre alignment works well
fit_adj = 0.5       #
fit_adjf = 0.1      #0.0
pre_range = 120     #points
fit_left = 70#20#38       #points
fit_right = 600#300     #points

Qr_0 = eval(sim_raw_head[9])
Qr_1 = eval(sim_raw_head[10])
Qr_p = sim_raw.shape[0]
Qr_array = np.linspace(Qr_0, Qr_1, Qr_p)


#----- Measured trace ---- (no loop required)
meas_data = np.array(meas_raw[1])
meas_data = norm_line(meas_data) #normalize
meas_time = np.array(meas_raw[0])#time in sec
meas_time = meas_time
#find zero position and adjust time
meas_pos = find_idx2(meas_data, fit_adjf, fit_adj, pre_range)
meas_time_off = meas_time[meas_pos]
meas_time = (meas_time - meas_time_off)
#crop data and time arrays
meas_data = crop_at_target(meas_data, meas_pos, fit_left, fit_right)
meas_time = crop_at_target(meas_time, meas_pos, fit_left, fit_right)
#normalize data for fitting
meas_data = norm_line(meas_data)


pl.figure(1)
pl.plot(meas_time,-meas_data+1) #plot trace to be fitted to

#----- Simulated trace ---- (loop is required)
sim_res_index = round(sim_raw.shape[1]/2) #selected freq index number to use
#produce time axis for simulated data set
sim_t_0 = eval(sim_raw_head[3]) #start in usec
sim_t_1 = eval(sim_raw_head[4])
sim_time2 = np.linspace(sim_t_0, sim_t_1, sim_raw.shape[2])
#create empty Xi square matrix
Ki2 = np.zeros(sim_raw.shape[0])
sim_store = np.zeros([sim_raw.shape[0],(fit_left+fit_right)])

for Qr_index in range(0,sim_raw.shape[0]):
    #Qr_index = 0 #for debugging
    sim_data = sim_raw[Qr_index][sim_res_index] #select Qr trace
    sim_data = sim_data[t_cut:]              #crop away edge
    sim_time = sim_time2[t_cut:]             #crop away edge
    #find zero position and adjust time
    sim_pos = find_idx2(sim_data, fit_adjf, fit_adj, pre_range)
    sim_time_off = sim_time[sim_pos]
    sim_time = (sim_time - sim_time_off)
    #crop data and time arrays
    sim_data = crop_at_target(sim_data, sim_pos, fit_left, fit_right)
    sim_time = crop_at_target(sim_time, sim_pos, fit_left, fit_right)
    #post normalize daya
    sim_data = norm_line(sim_data)
    sim_store[Qr_index] = sim_data 
    #calculate difference
    Ki2_tmp= (meas_data-sim_data)**2
    Ki2[Qr_index] = Ki2_tmp.sum()
    #plot results
    pl.figure(1)
    pl.plot(sim_time,-sim_data+1)


#plot Chi-squared as a func of Qr
pl.figure(2)
pl.plot(Qr_array, Ki2)
#plot best fit and measured data
pl.figure(3)
#pl.plot(meas_time,-meas_data+1)
pl.plot(-meas_data+1)
bestfit_idx = Ki2.argmin()
bestfit = -sim_store[bestfit_idx]+1
nodephasing = -sim_store[0]+1
#pl.plot(sim_time,bestfit)
pl.plot(bestfit)
#plot the difference squared (meas-bestfit)
pl.figure(4)
difference = ((-meas_data+1)-bestfit)
pl.plot(sim_time,difference**2)
print (difference**2).sum()



#save for gnuplot
#gnu_data1 = np.array(([meas_time,(-meas_data+1),bestfit, nodephasing]))
gnu_data1 = np.array(([meas_time,(-meas_data+1),bestfit])) # time, meas, fit
gnu_data2 = np.array(([Qr_array, Ki2]))

parser.savedat(filename_4, gnu_data1, delimiter = '\t')
parser.savedat(filename_5, gnu_data2, delimiter = '\t')