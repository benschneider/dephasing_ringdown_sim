# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 13:25:52 2014

@author: benschneider
just placing all parameters in here:
"""
import numpy as np
#import tkFileDialog
#file_path= tkFileDialog.askopenfilename()
#Nah!

#filename 1     input/output sim data file
#filename 2     input measured ringdown trace
#filename 3     input lockin filter response
#filename 4     output Best fitted trace
#filename 5     output Qr vs. Xi^2 values
filename_1 = 'sim_07850_6_im.mtx'
#'measured/07850_6_im_analysis_sub4.l.181.linecut.dat'
#'measured/07850_6_im_analysis_sub3.l.183.linecut.dat'
#'measured/07850_6_im_analysis_sub0.l.185.linecut.dat'
filename_2 = 'measured/07850_6_im_analysis_sub3.l.183.linecut.dat'
filename_3 = 'measured/07850_6_im_analysis_filt1.l.12.linecut.dat'
filename_4 = 'results/fit_07850_6_im.txt'
filename_5 = 'results/chi_07850_6_im.txt'

#filename_1 = 'sim_08267_8_im.mtx'
#filename_2 = 'measured/08267_8_im_analysis_sub5.l.130.linecut.dat'
#filename_3 = 'measured/08267_8_im_analysis_filt5.l.349.linecut.dat'
#filename_4 = 'results/fit_08267_8_im.txt'
#filename_5 = 'results/chi_08267_8_im.txt'

#filename_1 = 'sim_07225_3_im.mtx' 
#filename_2 = 'measured/07225_3_im_analysis_sub.l.260.linecut.dat' 
#filename_3 = 'measured/07225_3_im.l.0.linecut.dat' 
#filename_4 = 'results/fit_07225_3_im.txt' 
#filename_5 = 'results/chi_07225_3_im.txt'

#filename_1 = 'sim_08475_9_im.mtx' #sim data
#filename_2 = 'measured/08475_9_im_analysis_sub.l.109.linecut.dat'
#filename_3 = 'measured/08475_9_im_analysis.l.0.linecut.dat'
#filename_4 = 'results/fit_08475_9_im.txt'
#filename_5 = 'results/chi_08475_9_im_.txt'

''' this is for creating the simulated matrix '''
Qs = 3905           #Found Qs value
Qr_0 = Qs           #start spectral Q (excluding dephasing) 
Qr_1 = Qs+4000      #stop spectral Q (including dephasing)
Qr_p = 9            #num. of points
#Ringdowntime e.g. in usec (ensure same number of points per time!)
t_0 = 0             #start
t_1 = 71.7718       #stop
t_p = 1196          #num. of points
t_padd = 600        #add some extra points (uses same p/t) to avoid edge effects
t_cut = 50
t_p = t_p + t_padd
t_0 = t_0 - (t_1-t_0)/t_p*t_padd
#Frequency e.g. in MHz(add some extra points to avoid edge effects)
w_res = 302.545*2*np.pi #resonance position
span = 0.4*2*np.pi      #span around w_res/Qs*10
w_p = 1501              #number of points
w_0 = w_res-span        #start frequency (x2pi)
w_1 = w_res+span        #stop frequency (x2pi)