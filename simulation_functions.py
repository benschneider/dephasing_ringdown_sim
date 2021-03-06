# -*- coding: utf-8 -*-
"""
Created on Tue Apr 29 15:42:37 2014

@author: benschneider

- creates a 3d matrix.
"""
import scipy.ndimage
import scipy.signal as signal
from scipy.optimize import curve_fit
import equations as eq
import numpy as np


def get_ringdown_X(w_array, w_res, t, Q):
    X = np.zeros([len(w_array), len(t)]) #Get empty 2d Matrix
    ''' creates first matrix without any dephasing
    This creates the X quadrature output of a lockin measuring 
    the mechancal ringdown.
    '''
    i = 0
    for f in w_array:
        Xrd = eq.Xr(t, f, w_res, Q) #a ringdown list 10,9,8,7,5,4,2...0
        
        val1 = eq.Xr(0, f, w_res, Q) #eq.Yr(0, f, w_res, Q)
        j = 0
        while t[j] < 0:
            Xrd[j] = val1
            j += 1
              
        #Yrd = signal.convolve(Yrd,data_filt3) #convolve signal with lockin time response function
        X[i] = Xrd[0:len(t)] #assing to a 2d data map. (X vs Y)
        i += 1
    return X

def get_ringdown_X2(w_array, w_res, t, Qr):
    matrix3d = np.zeros((len(Qr), len(w_array), len(t))) #The result matrix    
    k = 0    
    for Q in Qr:
        matrix3d[k]  = get_ringdown_X(w_array, w_res, t, Q)
        k +=1
    return matrix3d

def get_ringdown_Y(w_array, w_res, t, Q, norm = False):
    Y = np.zeros([len(w_array), len(t)]) #Get empty 2d Matrix
    ''' creates first matrix without any dephasing
    This creates the Y quadrature output of a lockin measuring 
    the mechancal ringdown.
    '''
    i = 0
    for f in w_array:
        Yrd = eq.Yr(t, f, w_res, Q) #a ringdown list 10,9,8,7,5,4,2...0

        val1 = eq.Yr(0, f, w_res, Q) #eq.Yr(0, f, w_res, Q)
        j = 0
        while t[j] < 0:
            Yrd[j] = val1
            j += 1

        #Yrd = signal.convolve(Yrd,data_filt3) #convolve signal with lockin time response function
        Y[i] = Yrd[0:len(t)] #assing to a 2d data map. (X vs Y)
        i += 1

        if norm:
            Y = Y/Q #normalize ouput, make highest value = 1

    return Y

def get_ringdown_matrix_Y(w_array, w_res, t, Qr_array):
    matrix3d = np.zeros((len(Qr_array), len(w_array), len(t))) #The result matrix    
    k = 0
    for Q in Qr_array:
        matrix3d[k]  = get_ringdown_Y(w_array, w_res, t, Q)
        k +=1
    return matrix3d
#-------- Convolution function
#-------- calculate dephasing     


def get_dephased_matrix(matrix3d, Qd_array, w_array, w_res, method='gaus'):
    ''' matrix is a 3d np.array, have 2-d maps which require different dephasing 
        method can be gausian dephasing or resonator lorenzian like dephasing.
        method = 'gaus' > gaussian dehasing
        method = 'lor'  > lor like shape dephasing
    '''
    if method == 'lor':
        matrix3d = _get_dephased_matrix_lor(matrix3d, Qd_array, w_res, w_array)
    elif method == 'gaus':
        matrix3d = _get_dephased_matrix_gaus(matrix3d, Qd_array, w_res, w_array)
    return matrix3d
    
def _get_dephased_matrix_lor(matrix3d, Qd_array, w_res, w_array): 
    ''' dephasing convolved with a resonator like line shape'''
    k = 1
    w_res = w_array.mean() #use center position for dephasing calculation
    for Qd in Qd_array:
        filter2d =  filter2d_lor(w_array, w_res, Qd)#,crop = int(matrix3d.shape[1]/2.05))    
        matrix3d[k] = signal.convolve2d(matrix3d[k], filter2d, mode = 'same') #and here python commits suicide...
        #matrix3d[k] = matrix3d[k]/adj #for debuging of the convolution function 
        k += 1
        print str(100.0*k/(matrix3d.shape[0])) + ' %'
    return matrix3d

def _get_dephased_matrix_gaus(matrix3d, Qd_array, w_res, w_array):
    ''' dephasing convolved with a gaussian line shape
    This function is a lot faster!    
    '''
    k = 1
    #w_res = w_array.mean() #resonance frequency    
    fpix = eq.pix2f(w_array) #(pixel/frequency)
    for Qd in Qd_array:
        FWHM = w_res/Qd #width in freq for given Q
        pixel = FWHM*fpix/2  #in terms of number of pixel
        matrix3d[k] = scipy.ndimage.filters.gaussian_filter(matrix3d[0], (pixel, 0))
        k += 1
    return matrix3d

def get_matrix_lockin_convolved(matrix3d, filt_fun):
    '''convolve a 3d-matrix with a function'''
    for k in range(0,matrix3d.shape[0]):
        matrix3d[k] = convolution_1d(matrix3d[k], filt_fun)
        #signal.convolve2d(matrix3d[k], filtfun, mode = 'same') #AAAARGGHHH!!!!.... this takes too long
        print str(100.0*k/(matrix3d.shape[0]-1)) + ' %'
    return matrix3d


#-------- Filter functions

def filter2d_lockin(lockin_response):
    '''
    is used to extract lockin-filter
    '''
    #norm signal
    lockin_output = lockin_response - lockin_response.min()
    lockin_output = lockin_output/lockin_output.max()
    raw_filter = diff_acurate(lockin_output) #calc derivative
    return norm_filter(raw_filter[1:-1])
    


#    tmp = np.zeros([1,len(raw_filter)]) #store filter in the right format
#    for j in range(0, tmp.shape[0]):
#        tmp[j] = norm_filter(raw_filter) #store normalized filter


def diff_acurate(a):
    '''
    create an position accurate forward derivative (using 3 points)    
    this works on 1 and 2 dim arrays
    takes the derivative from left to right
    '''
    tmp_m = np.concatenate( ([ [a.transpose()[0]], a.transpose(), [a.transpose()[-1]]  ]) ,axis = 0)
    tmp_1 = np.concatenate( ([ [a.transpose()[0]], [a.transpose()[0]], a.transpose()   ]) ,axis = 0)
    tmp_2 = np.concatenate( ([ a.transpose(), [a.transpose()[-1]], [a.transpose()[-1]] ]) ,axis = 0)
    tmp_diff = ((tmp_2 - tmp_m) + (tmp_m -tmp_1))/2.0
    return tmp_diff.transpose()


# for a range of Qr and a fixed assumed Qs
def filter2d_lor(w_array, w_res, Qd, crop = 0):
    ''' -
    w_array is the list of frequencies
    Qd the dephasing Q factor 
    if crop > 0 it crops from both ends points by numbers given
    '''
    
    function = eq.yqfit(w_array, Qd, w_res) #obtain a lorenzian line shape
    function = -function/Qd
    function = norm_filter(function) #normalize filter
    #this simply crops off the edges of the filter    
    if crop > 0:
        function = function[crop:-crop]
    
    filter2d = np.zeros([1,len(function)])

    for j in range(0, filter2d.shape[0]):
        filter2d[j] = function
        
    filter2d = filter2d.transpose()
    return filter2d


def norm_filter(filterfun):
    '''
    go across each element and adds them
    this value is then subtracted from the filter function
    (N = integral{(-inf) to (+inf)}{function})
    '''
    N = 0
    for i in range(0,len(filterfun)):
        N += filterfun[i]

    return filterfun/N


#-------- for testing
def get_normed_matrix(matrix3d, Qr):
    '''
    this function simply normalized the matrix using its Q-factors
    '''
    k = 0    
    for Q in Qr:
        matrix3d[k] = matrix3d[k] / Q
        k += 1
    return matrix3d
#------- Convolution functions

def convolution_1d(data2d, filt_fun):
    tmp = data2d
    for j in range(0,len(tmp)):
        tmp[j] = _convolution(data2d[j], filt_fun)
    return tmp

def _convolution(data, filt_func):
    '''
    This function does calculate a 1D convolution
    ,the goal is to test a convolution method.
    (it is not yet optimized for speed)
    expected input: func is the data to be processed 
    and filt_func is the filter to be used
    '''
    length = len(data)
    filt_new = _get_convolution_filter(filt_func,length)
    #new_filt =  [0,0,0,0,0,0,0.5,1,0.5,0,0,0,0,0,0]
    filt_length = len(filt_func)
    
    #dim1 = (len(data)+filt_length)
    
    dim2 = len(data)
    test_array = np.zeros([dim2,dim2])

    for i in range(0,dim2):
        #for each element in data do
        filt_start = length+filt_length/2-i
        filt_stop  = (2*length+filt_length/2-i)
        filt_window = filt_new[filt_start:filt_stop]
        #print filt_window        
        test_array[i] = filt_window*data[i]
    
    data_convolved = test_array.sum(0) 
    return data_convolved
    
def _get_convolution_filter(filt_func,length):
    '''puts the filter in the center and adds sufficient zeros at either side'''
    #length = len(filt_func)
    a = np.array([0]*length)
    b = np.array([0]*(length+1))    
    new_filt = np.concatenate((a,filt_func,b), axis = 0)
    return new_filt

#---------
def fit_mat_ringdown(t, matrix3d):
    '''  Fit the created ringdown matrix, 
    it will assume a matrix with the following shape and form: 
    and then return a result matrix
    '''
    num_dephasing = matrix3d.shape[0]
    num_freq = matrix3d.shape[1]
    
    #get position where t = 0
    t0_index = int(-t[0]/((t[-1] - t[0])/matrix3d.shape[2])) 

    rpopt = np.zeros([num_dephasing, 3])
    rpcov = np.zeros([num_dephasing, 3, 3])
    qrfit = np.zeros([1, num_dephasing*2, len(t)])

    for k in range(0, num_dephasing):
        #Fit data to eq.expfit
        Ringdown = matrix3d[k][num_freq/2]
        rpopt[k], rpcov[k] = curve_fit(eq.expfit, t[t0_index:-1], Ringdown[t0_index:-1],
                    p0=(6000, -6000, 1), sigma=None, maxfev=5000)

        # create fitted function
        fited = eq.expfit(t, *rpopt[k])
        fited[0:t0_index] = fited[t0_index]
        
        #store data plus fit lines into a matrix
        qrfit[0][k*2] = Ringdown
        qrfit[0][k*2+1] = fited

    Qrs = rpopt.transpose()
    return Qrs, qrfit
   
def fit_matrix_spectral(w_array, matrix3d):
    '''  Fit the created ringdown matrix, 
    it will assume a matrix with the following shape and form: ...
    output:
    fit results Q factor per dephasing
    matrix containing the fit and dephasing positions used
    '''
    num_dephasing = matrix3d.shape[0]
    spopt = np.zeros([num_dephasing, 1])
    spcov = np.zeros([num_dephasing, 1, 1])
    qsfit = np.zeros([1, num_dephasing*2, len(w_array)])

    for k in range(0, num_dephasing):
        #Fit data to eq.yfit
        temp3 = matrix3d[k].transpose()
        Spectral = temp3[0]
        spopt[k], spcov[k] = curve_fit(eq.yqfit, w_array, Spectral, p0=1500, 
                                                    sigma=None, maxfev=5000)
        
        #Create fitted function:
        fited = eq.yqfit(w_array, spopt[k])
        
        #store data plus fit lines
        qsfit[0][k*2] = Spectral
        qsfit[0][k*2+1] = fited

    Qsp = spopt.transpose()   
    return Qsp, qsfit
    
def shrink_extend_array(a, xdim):
    ''' 
    this function is used to shrink or extend an 1d-array
    as a first approximation, it takes the average between 2 points! 
    to extend or to shrink the array
    so any changes larger than a factor of 2 or 3 will be inaccurate
    '''
    b = np.zeros((xdim))
    
    
    for i in range(xdim):
        index = (i * len(a) / xdim)
        #the would be index     

        tmp1 = a[int(index)] 
        tmp2 = a[round(index)]
        tmp = (tmp1+tmp2)/2
        #tmp = a[index]
        
        b[i] = tmp
    
    return b