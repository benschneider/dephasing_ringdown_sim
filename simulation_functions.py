# -*- coding: utf-8 -*-
"""
Created on Tue Apr 29 15:42:37 2014

@author: benschneider

- creates a 3d matrix.
"""
import scipy.ndimage #a fast filter
from scipy.optimize import curve_fit
import equations as eq
import numpy as np

def get_ringdown_X(w, w0, t, Q):
    X = np.zeros([len(w), len(t)]) #Get empty 2d Matrix
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
              
        #Yrd = signal.convolve(Yrd,data_filt3) #convolve signal with lockin time response function
        X[i] = Xrd[0:len(t)] #assing to a 2d data map. (X vs Y)
        i += 1
    return X
    
def get_ringdown_Y(w, w0, t, Q):
    Y = np.zeros([len(w), len(t)]) #Get empty 2d Matrix
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
        #Yrd = signal.convolve(Yrd,data_filt3) #convolve signal with lockin time response function
        Y[i] = Yrd[0:len(t)] #assing to a 2d data map. (X vs Y)
        i += 1
    return Y

def get_dephased_matrix(matrix, dephasing):
    k = 1
    matrix.shape[1]
    matrix.shape[2]
    for d in dephasing[1:]: #d = sigma
        #tmp_mat = np.zeros((len(t), len(w)))
        #l = 0
        #do a first order LP req 2usec see paper by Young & Vliet
        matrix[k] = scipy.ndimage.filters.gaussian_filter(matrix[0], (d, 0))
        k += 1
    return matrix

def fit_mat_ringdown(t, matrix):
    '''  Fit the created ringdown matrix, 
    it will assume a matrix with the following shape and form: 
    and then return a result matrix
    '''
    num_dephasing = matrix.shape[0]
    num_freq = matrix.shape[1]

    rpopt = np.zeros([num_dephasing, 3])
    rpcov = np.zeros([num_dephasing, 3, 3])
    qrfit = np.zeros([1, num_dephasing*2, len(t)])

    for k in range(0, num_dephasing):
        #Fit data to eq.expfit
        Ringdown = matrix[k][num_freq/2]
        rpopt[k], rpcov[k] = curve_fit(eq.expfit, t[0:400], Ringdown[0:400],
                    p0=(6000, -6000, 1), sigma=None, maxfev=5000)
        #pl.pcolor(pcov)
    
        #store data plus fit lines
        qrfit[0][k*2] = Ringdown
        qrfit[0][k*2+1] = eq.expfit(t, *rpopt[k])

    Qrs = zip(*rpopt)
    return Qrs, qrfit
   
def fit_matrix_spectral(w, matrix):
    '''  Fit the created ringdown matrix, 
    it will assume a matrix with the following shape and form: 
    output:
    fit results Q factor per dephasing
    matrix containing the fit and dephasing positions used
    '''
    num_dephasing = matrix.shape[0]
    spopt = np.zeros([num_dephasing, 1])
    spcov = np.zeros([num_dephasing, 1, 1])
    qsfit = np.zeros([1, num_dephasing*2, len(w)])

    for k in range(0, num_dephasing):
        #Fit data to eq.yfit
        temp3 = zip(*matrix[k])
        Spectral = temp3[0]
        spopt[k], spcov[k] = curve_fit(eq.yqfit, w, Spectral, p0=6000, 
                                                    sigma=None, maxfev=5000)
        
        #store data plus fit lines
        qsfit[0][k*2] = Spectral
        qsfit[0][k*2+1] = eq.yqfit(w, spopt[k])

    Qsp = zip(*spopt)    
    return Qsp, qsfit

