# -*- coding: utf-8 -*-
"""
Created on Tue Apr 29 15:42:37 2014

@author: benschneider

- creates a 3d matrix.
"""

import equations as eq
import numpy as np

def get_ringdown_X(w, t, f, w0, Q):
    X = np.zeros([len(w), len(t)]) #Yes lets chew up some temp.memory
    ''' creates first matrix without any dephasing
    This creates the X quadrature output of a lockin measuring 
    the mechancal ringdown.
    '''
    i = 0
    for f in w:
        Xrd = eq.Xr(t, f, w0, Q) #a ringdown list 10,9,8,7,5,4,2...0
        
        val1 = eq.Xr2(0, f, w0, Q) #eq.Yr(0, f, w0, Q)
        j = 0
        while t[j] < 0:
            Xrd[j] = val1
            j += 1
              
        #Yrd = signal.convolve(Yrd,data_filt3) #this can be uncommented to deactivate the lockin filter sim.
        X[i] = Xrd[0:len(t)] #assing to a 2d data map. (X vs Y)
        i += 1
    return X
    
def get_ringdown_Y(w, t, f, w0, Q):
    Y = np.zeros([len(w), len(t)]) #Yes lets chew up some temp.memory
    ''' creates first matrix without any dephasing
    This creates the Y quadrature output of a lockin measuring 
    the mechancal ringdown.
    '''
    i = 0
    for f in w:
        Yrd = eq.Yr(t, f, w0, Q) #a ringdown list 10,9,8,7,5,4,2...0
        
        val1 = eq.Yr2(0, f, w0, Q) #eq.Yr(0, f, w0, Q)
        j = 0
        while t[j] < 0:
            Yrd[j] = val1
            j += 1
              
        #Yrd = signal.convolve(Yrd,data_filt3) #this can be uncommented to deactivate the lockin filter sim.
        Y[i] = Yrd[0:len(t)] #assing to a 2d data map. (X vs Y)
        i += 1
    return Y
