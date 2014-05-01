# -*- coding: utf-8 -*-
"""
Created on Tue Apr 29 15:42:37 2014

@author: benschneider

- creates a 3d matrix.
"""
import scipy.ndimage #a fast filter
from scipy.optimize import curve_fit
import scipy.signal as signal
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

def get_ringdown_X2(w, w0, t, Qr):
    matrix = np.zeros((len(Qr), len(w), len(t))) #The result matrix    
    k = 0    
    for Q in Qr:
        matrix[k]  = get_ringdown_X(w, w0, t, Q)
        k +=1
    return matrix

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

def get_ringdown_Y2(w, w0, t, Qr):
    matrix = np.zeros((len(Qr), len(w), len(t))) #The result matrix    
    k = 0
    for Q in Qr:
        matrix[k]  = get_ringdown_Y(w, w0, t, Q)
        k +=1
    return matrix

# OLD NOT IN USE
#def get_dephased_matrix(matrix, dephasing):
#    k = 1
#    for d in dephasing: #d = sigma
#        #do a first order LP req 2usec see paper by Young & Vliet
#        matrix[k] = scipy.ndimage.filters.gaussian_filter(matrix[0], (d, 0))
#        k += 1
#    return matrix

def get_dephased_matrix(matrix, dephasing, w, method='gaus'):
    ''' method can be gausian dephasing or resonator lorenzian like dephasing.
        method = 'gaus' > gaussian dehasing
        method = 'lor'  > lor like shape dephasing
    '''
    if method == 'lor':
        matrix = _get_dephased_matrix_lor(matrix, dephasing, w)
    elif method == 'gaus':
        matrix = _get_dephased_matrix_gaus(matrix, dephasing, w)
    return matrix
    
def _get_dephased_matrix_lor(matrix, dephasing, w):    
    ''' dephasing convolved with a resonator like line shape'''
    k = 1
    for Qd in dephasing:
        filter2d =  get_filter2d(w,Qd, crop = int(matrix.shape[1]/2.1))    
        matrix[k] = signal.convolve2d(matrix[k], filter2d, mode = 'same') #and here python commits suicide...
        k += 1
    return matrix

def _get_dephased_matrix_gaus(matrix, dephasing, w):
    ''' dephasing convolved with a gaussian line shape
    This function is a lot faster!    
    '''
    k = 1
    w0 = w.mean() #resonance frequency    
    fpix = eq.pix2f(w) #(pixel/frequency)
    for Qd in dephasing:
        FWHM = w0/Qd #width in freq for given Q
        pixel = FWHM*fpix/2  #in terms of number of pixel
        matrix[k] = scipy.ndimage.filters.gaussian_filter(matrix[0], (pixel, 0))
        k += 1
    return matrix

# for a range of Qr and a fixed assumed Qs
def get_filter2d(w, Qd, crop = 0):
    ''' -
    w is the list of frequencies
    Qd the dephasing Q factor 
    and crop != 0 crops from both ends points
    '''
    function = eq.yqfit(w, Qd) #obtain a lorenzian line shape
    function = -function/Qd
    #this simply crops off the edges of the filter    
    if crop > 0:
        function = function[crop:-crop]
    
    filter2d = np.zeros([1,len(function)])

    for j in range(0, filter2d.shape[0]):
        filter2d[j] = function
        
    filter2d = zip(*filter2d) #rotate
    return filter2d


def fit_mat_ringdown(t, matrix):
    '''  Fit the created ringdown matrix, 
    it will assume a matrix with the following shape and form: 
    and then return a result matrix
    '''
    num_dephasing = matrix.shape[0]
    num_freq = matrix.shape[1]
    
    #get position where t = 0
    t0_index = int(-t[0]/((t[-1] - t[0])/matrix.shape[2])) 

    rpopt = np.zeros([num_dephasing, 3])
    rpcov = np.zeros([num_dephasing, 3, 3])
    qrfit = np.zeros([1, num_dephasing*2, len(t)])

    for k in range(0, num_dephasing):
        #Fit data to eq.expfit
        Ringdown = matrix[k][num_freq/2]
        rpopt[k], rpcov[k] = curve_fit(eq.expfit, t[t0_index:-1], Ringdown[t0_index:-1],
                    p0=(6000, -6000, 1), sigma=None, maxfev=5000)

        # create fitted function
        fited = eq.expfit(t, *rpopt[k])
        fited[0:t0_index] = fited[t0_index]
        
        #store data plus fit lines into a matrix
        qrfit[0][k*2] = Ringdown
        qrfit[0][k*2+1] = fited

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
        
        #Create fitted function:
        fited = eq.yqfit(w, spopt[k])
        
        #store data plus fit lines
        qsfit[0][k*2] = Spectral
        qsfit[0][k*2+1] = fited

    Qsp = zip(*spopt)    
    return Qsp, qsfit

