'''
Ben Schneider

this script generates a simulated ringdown dataset
be carefull not to increase the limits too much a high
(in consumes space and memory)

import numpy as np
'''
#from numpy import sin, cos, e, pi, sqrt, linspace, array, zeros, ones, savetxt
import numpy as np
#from time import time
from scipy.optimize import curve_fit
from scipy import signal
#import scipy.signal as signal
#from scipy.signal.convolve
import scipy.ndimage #The fastest filter
import pylab as pl

#import matplotlib.pyplot as pl


#t0 = time()
#execfile('mtx_parser.py') #to save data as mtx file
#execfile('gnu.py') # to plot using gnuplot
#execfile('equations.py') # load my equations
import equation as eq
import mtxparser as mp



filename = 'sim_dephasing_lockin.mtx'
Q = 6000.0
w0 = 300.0*2*np.pi
dephasing = np.linspace(0, 50, 50)#in (Mhz*2.355) (FWHM)

span = w0/Q*30 #0.1*pi;
rdtime = int(Q/w0*12.0)
t = np.linspace(0.0, rdtime, 601) #Ringdowntime in usec
w = np.linspace(w0-span, w0+span, 601) #Frequency range

x = eq.Xd(w, w0, Q)
y = eq.Yd(w, w0, Q)

pl.plot(w/2/np.pi, x)
pl.plot(w/2/np.pi, y)
#pl.show()

#gn.plot(Gnuplot.Data(w/2/pi,x))
#gn.replot(Gnuplot.Data(w/2/pi,y))

Y = np.zeros([len(w), len(t)]) #Yes lets chew up some temp.memory
matrix = np.zeros((len(dephasing), len(w), len(t))) #The result matrix

#For ringdown fits
rpopt = np.zeros([len(dephasing), 3])
rpcov = np.zeros([len(dephasing), 3, 3])
spopt = np.zeros([len(dephasing), 1])
spcov = np.zeros([len(dephasing), 1, 1])

#Matrix containing the linecut and its fit to be seen in spyview
qsfit = np.zeros([1, len(dephasing)*2, len(w)])
qrfit = np.zeros([1, len(dephasing)*2, len(t)])

#extract the lockin convolution filter from (a step function)
mtx1_data, mtx1_header = mp.loadmtx('meas_data/n20dBm_amp.mtx')
mtx1_data = mtx1_data[0] #it is a 2d colour map no 3d req.
data_off = mtx1_data[0] - mtx1_data[0].min()
data_off = data_off/data_off.max()
data_filt = -np.diff(data_off,1) #obtain Filter
data_filt2 = data_filt[2000:4000] #croped for this particular data set

#some experimental averaging
fact = 1
a = np.zeros(len(t)/fact)
a = np.array(a)
avg = float(len(data_filt))/float(len(t))*fact
avg = int(avg)
j = 0
j2 = 0
for i in range(0,len(data_filt)):
    j = i%avg
    if j == 0 and i != 0:
        j2 +=1
    if (j+j2) < len(t)/fact:
        a[j+j2] += data_filt[i]/avg
    else:
        break;
data_filt3 = a[170:330] #crop an averaged filter

#creates first matrix without any dephasing
i = 0
for f in w:
    Yrd = eq.Yr(t, f, w0, Q) #a ringdown list 10,9,8,7,5,4,2...0
    Yrd = signal.convolve(Yrd,data_filt3) #this can be uncommented to deactivate the lockin filter sim.
    Y[i] = Yrd[0:len(t)] #assing to a 2d data map. (X vs Y)
    i += 1
matrix[0] = Y #asign to a 3d data map (X vs Y vs Z)

#do lowpass of the matrix
k = 1
for d in dephasing[1:]: #d = sigma
    temp2 = np.zeros((len(t), len(w)))
    l = 0
    #do a first order LP req 2usec see paper by Young & Vliet
    matrix[k] = scipy.ndimage.filters.gaussian_filter(matrix[0], (d, 0))
    k += 1

#fit the first stuff
#optimize.leastsq
#temp3 = zip(*matrix[0])
#Spectral = temp3[0]
#spopt[0], spcov[0] = curve_fit(yqfit, w, Spectral,
#                    p0=6000, sigma=None, maxfev=5000)
#
#Ringdown = matrix[0][len(w)/2]
#rpopt[0], rpcov[0] = curve_fit(expfit, t, Ringdown,
#                    p0=(6000, -6000, 0), sigma=None, maxfev=5000)
#
##store data plus fit lines
#qsfit[0] = Spectral
#qsfit[1] = yqfit(w,spopt[0])
#qrfit[0] = Ringdown
#qrfit[1] = expfit(t,*rpopt[0])

#fit dephased data.
for k in range(0, len(dephasing)):
    #Fit data:
    temp3 = zip(*matrix[k])
    Spectral = temp3[0]
    spopt[k], spcov[k] = curve_fit(eq.yqfit, w, Spectral,
                p0=6000, sigma=None, maxfev=5000)

    Ringdown = matrix[k][len(w)/2]
    rpopt[k], rpcov[k] = curve_fit(eq.expfit, t[0:400], Ringdown[0:400],
                p0=(6000, -6000, 1), sigma=None, maxfev=5000)
    #pl.pcolor(pcov)

    #store data plus fit lines
    qsfit[0][k*2] = Spectral
    qsfit[0][k*2+1] = eq.yqfit(w, spopt[k])
    qrfit[0][k*2] = Ringdown
    qrfit[0][k*2+1] = eq.expfit(t, *rpopt[k])


    #plot results
    pl.subplot(3, 1, 1)
    pl.plot(w/2/np.pi, Spectral)
    pl.plot(w/2/np.pi, eq.yqfit(w, spopt[k][0]))
    pl.hold(True)
    pl.subplot(3, 1, 2)
    pl.plot(t, Ringdown)
    pl.plot(t, eq.expfit(t, rpopt[k][0], rpopt[k][1], rpopt[k][2]))
    pl.hold(True)

Qrs = zip(*rpopt)
Qsp = zip(*spopt)
pl.subplot(3, 1, 3)
pl.plot(dephasing, Qrs[0])
pl.plot(dephasing, Qsp[0])
pl.hold(True)


 #Dephasing (per pixel) per (w0/Q)
deph_fact = ((w.max()-w.min())/len(w))*Q/w0

#prepare stuff to be saved into the txt file
stuff = (dephasing*deph_fact, Qrs[0], Qsp[0])
mp.savedat("QsQr.txt", stuff, delimiter='\t')

#Save the Spectral fitting matix into one MTX file
head = ['Units', 'Qs_fits',
        'Dephasing (rel to FWHM)', '0', str(deph_fact*dephasing[-1]),
        'RF frequency (MHz)', str(w[-1]/2/np.pi), str(w[0]/2/np.pi),
        'other', '0', '1']

mp.savemtx('qsfit.mtx', eq.qsfit, header=head) #save as MTX

#Save the Ringdown fitting matix into one MTX file
head = ['Units', 'Qs_fits',
        'Dephasing (rel to FWHM)', '0', str(deph_fact*dephasing[-1]),
        'Time (us)', str(t[-1]), str(t[0]),
        'other', '0', '1']
mp.savemtx('qrfit.mtx', qrfit, header=head) #save as MTX

#Save the dephasing matix into another MTX file
head = ['Units', 'X Y R [V]#values',
        'Time (us)', str(t[0]), str(t[-1]),
        'RF frequency (MHz)', str(w[-1]/2/np.pi), str(w[0]/2/np.pi),
        'Dephasing (rel to FWHM)', '0', str(deph_fact*dephasing[-1])]
mp.savemtx(filename, matrix, header = head) #save as MTX
