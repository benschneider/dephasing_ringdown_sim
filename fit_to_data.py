# -*- coding: utf-8 -*-
"""
Created on Tue May  6 17:54:15 2014

@author: benschneider

structure of this code:
1st.
load simulated data with the fixed Qs and different Qrs.
load measured data with known Qs and unknown Qr 
(and pick the right linecut)

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

#filenames
filename1 = 'simulated data set'
filename2 = 'measured data set'

#load data
data_sim,  data_sim_head  = parser.loadmtx(filename1) 
data_meas, data_meas_head = parser.loadmtx(filename2)


# part of picking Qrs
ringdown_meas = data_meas[100] #pick the right on resonance line number

for Qr in range(0,len(data_sim)):
    ringdown_sim = data_sim[Qr]

