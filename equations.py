#from numpy import sin, cos, e, sqrt, pi
import numpy as np
#def hho(w,w0,q):
#  return w0*w0/(w0*w0-w*w + i*w*w0/Q)

def Xd(w,w0,Q):
  return (w0*w0*w0*w0-w0*w0*w*w)/( (w0*w0-w*w)*(w0*w0-w*w) + w*w*w0*w0/Q/Q )

def Yd(w,w0,Q):
  return -w0*w0*w0*w/Q/( (w0*w0-w*w)*(w0*w0-w*w) + w*w*w0*w0/Q/Q )

def Xr(t,w,w0,Q):
    return np.e**(-w*t/(2*Q))*(Xd(w,w0,Q)*np.cos(t*(w0 - w)) + Yd(w,w0,Q)*w/w0*np.sin(t*(w0 - w)))

def Yr(t,w,w0,Q):
    return np.e**(-w*t/(2*Q))*(-Xd(w,w0,Q)*np.sin(t*(w0 - w)) + Yd(w,w0,Q)*w/w0*np.cos(t*(w0 - w)))

def Xr2(t,w,w0,Q):
    return np.e**(-w*t/(2*Q))*(Xd(w,w0,Q)*np.cos(t*(w0 - w))**2 + Yd(w,w0,Q)*w/w0*np.sin(t*(w0 - w))**2)

def Yr2(t,w,w0,Q):
    return np.e**(-w*t/(2*Q))*(-Xd(w,w0,Q)*np.sin(t*(w0 - w))**2 + Yd(w,w0,Q)*w/w0*np.cos(t*(w0 - w))**2)

def yqfit(w, Q, w0 =300*2*np.pi):
    ''' Equation for fitting (Spectral Q)
    returns -w0*w0*w0*w/Q/( (w0*w0-w*w)*(w0*w0-w*w) + w*w*w0*w0/Q/Q )    
    '''
    #w0 = 300*2*np.pi; #at a fixed freq
    return -w0*w0*w0*w/Q/( (w0*w0-w*w)*(w0*w0-w*w) + w*w*w0*w0/Q/Q )

def expfit(t, Q, a, b, w0 =300*2*np.pi):
    '''Equation for fitting (Ringdown Q)'''
    #w0 = 300*2*np.pi; #at a fixed freq
    return (a*np.e**(-t*w0/(2*Q))+b)


def Qd(Qs,Qr):
    ''' Equation returns an estimate for Qr to be used, 
    given an estimate to Qr for a given Qs and Qd 
    '''
    tmp = 1.0/(1.0/float(Qs)-1.0/Qr)
    return tmp


def pix2f(w):
    '''
    this little handy fuction simply returns a factor such that one can easily 
    convert i.e. a frequency to number of points
    it requires the array of frequency
    then from the frequency span and number of points
    it returns a factor describing the number of points required to change frequency by #1    
    '''
    w_span = w[-1]-w[0] #frequency range *2pi
    w_p = len(w) #number of pixel
    fpix = w_p/w_span #pixel required to shift freq by 1 MHz
    return fpix


#Equation for lowpass (returns a gaussian)
#Doing the Lowpas in C,
#Gaussian filtering can be implemented using recursion (Young & Vliet,
#1995; Signal Processing, 44: 139-151). The recursive approximation is
#very accurate and does not introduce ringing. It is anti-causal
#(forward-backward) and has zero phase response.

#or using a correlation funcion from the scipy lib.
# which is a lot slower but good enough for now.

'''
def gaus(sigma=1,pos = 50, w = range(0,100)): 
    g1 = [1 / (sigma * np.sqrt(2*np.pi)) * np.e**(-float(x-pos)**2/(2*sigma**2)) for x in w]
    #g2 = np.ones([m,len(w)]);
    return g1 #zip(*g1*g2)
    #a gaussian 


# the following can be derived with a simple taylor expansion of a gaussian
# and doing some fancy stuff with it which is cool but ... nah...

from numpy import array, zeros, ones, flipud, fliplr
from scipy.signal import lfilter
from math import sqrt

def __gausscoeff(s):
    if s < .5: raise ValueError, \
    'Sigma for Gaussian filter must be >0.5 samples'
    q = 0.98711*s - 0.96330 if s > 0.5 else 3.97156 \
    - 4.14554*sqrt(1.0 - 0.26891*s)
    b = zeros(4)
    b[0] = 1.5785 + 2.44413*q + 1.4281*q**2 + 0.422205*q**3
    b[1] = 2.44413*q + 2.85619*q**2 + 1.26661*q**3
    b[2] = -(1.4281*q**2 + 1.26661*q**3)
    b[3] = 0.422205*q**3
    B = 1.0 - ((b[1] + b[2] + b[3])/b[0])
    # convert to a format compatible with lfilter's
    # difference equation
    B = array([B])
    A = ones(4)
    A[1:] = -b[1:]/b[0]
    return B,A
    
def Gaussian1D(signal, sigma, padding=0):
    n = signal.shape[0]
    tmp = zeros(n + padding)
    if tmp.shape[0] < 4: raise ValueError, \
    'Signal and padding too short'
    tmp[:n] = signal
    B,A = __gausscoeff(sigma)
    tmp = lfilter(B, A, tmp)
    tmp = tmp[::-1]
    tmp = lfilter(B, A, tmp)
    tmp = tmp[::-1]
    return tmp[:n]
    
def Gaussian2D(image, sigma, padding=0):
    n,m = image.shape[0],image.shape[1]
    tmp = zeros((n + padding, m + padding))
    if tmp.shape[0] < 4: raise ValueError, \
    'Image and padding too small'
    if tmp.shape[1] < 4: raise ValueError, \
    'Image and padding too small'
    B,A = __gausscoeff(sigma)
    tmp[:n,:m] = image
    tmp = lfilter(B, A, tmp, axis=0)
    tmp = flipud(tmp)
    tmp = lfilter(B, A, tmp, axis=0)
    tmp = flipud(tmp)
    tmp = lfilter(B, A, tmp, axis=1)
    tmp = fliplr(tmp)
    tmp = lfilter(B, A, tmp, axis=1)
    tmp = fliplr(tmp)
    return tmp[:n,:m]
'''