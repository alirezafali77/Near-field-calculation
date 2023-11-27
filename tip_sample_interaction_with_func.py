# -*- coding: utf-8 -*-
"""
Created on Sun Jan 22 17:47:19 2023

@author: alire
"""

import matplotlib.pyplot as plt
import numpy as np
import math
import cmath   #import cmath for complex number operations
from scipy.integrate import quad
import mpmath as mp
from mpmath import *
"""
This code is to simulate near-field interaction of tip and sample for 4 layers of materials
first layer is air which always use 1 as dielectric function
second layer defined as "data_top" 
third layer defined as "data_substrate"
"data_top" and "data_substrate" are inputs that user needs to provide.
For this simulation, forth layer defined as Silicon; dielectric=11.7 "eps_sub_sub=11.7"
Important note: model is based on dielctric function. if data is based on refractive index, it has to convert
Equations have been written in the following paper:
https://pubs.acs.org/doi/abs/10.1021/acsphotonics.0c00951    
"""
#Variables
thickness_top=100e-9 #thickness of material on top
thickness_sub=300e-9  #thickness of first substrate
thickness_sub_sub=1000e-9  #thickness of main substrate
eps_sub_sub=11.7  #dielectric of silicon
q_factor=1e7
tip_radius=10e-9
tip_height=15*tip_radius
tip_distance=1e-9
harmonic=2  
probe_fre=260e3
time_tap=2*math.pi/probe_fre
number_of_steps=1000
max_range_of_cal=1500
min_range_of_cal=700
step_define=np.linspace(min_range_of_cal,max_range_of_cal,number_of_steps)

# dielectric information of materials
data_substrate=np.loadtxt('C:/Users/alire/OneDrive - University of Georgia/Desktop/python practice/calculation simulation/refractive_sio2.txt')
data_top=np.loadtxt('C:/Users/alire/OneDrive - University of Georgia/Desktop/python practice/calculation simulation/refractive_Mica.txt')

size_substrate=np.shape(data_substrate)
size_top=np.shape(data_top)

"""
Section A
change the substrate refractive index to dielectric function, 
this part is only because the dielctric informatin of the third layer is based on 
wavelenght, and refractive index 
"""

x2=10000/data_substrate[:,0]  #convert wavelength to wavenumber
y2=data_substrate[:,1]
z2=data_substrate[:,2]
eps_real=(y2*y2)-(z2*z2)
eps_imag=2*y2*z2
eps_substrate1=eps_real+1j*eps_imag
min_val=min(x2)
max_val=max(x2)
print('Range of the substrate is: [', min_val,',',max_val,']')
"""
section A ended
"""


# organize the top material dielectric function
x1=data_top[:,0]
eps_top1=data_top[:,1]+1j*data_top[:,2]

#
min_val_top=min(data_top[:,0])
max_val_top=max(data_top[:,0])
print('Range of top material is: [' , min_val_top, ',', max_val_top, ']')


"""
Section B
reproduce dielectrics in equal step size, 
Since the experimental data for the two material have different steps,
the following steps will make the diecltric functions in simillar steps
"""

xx_top=np.ones(number_of_steps)
eps_top=xx_top+1j*xx_top
eps_sub=xx_top+1j*xx_top
for x_step in range (number_of_steps):
    difference_array = np.absolute(x1-step_define[x_step])
    index = difference_array.argmin()
    eps_top[x_step]=data_top[index,1]+1j*data_top[index,2]
    difference_array_sub = np.absolute(x2-step_define[x_step])
    index_sub = difference_array_sub.argmin()
    eps_sub[x_step]=eps_real[index_sub]+1j*eps_imag[index_sub]
    
"""
section B ended
"""   

# calculations for the model
g=.7*cmath.exp(1j*.06)
C=tip_height*tip_radius*tip_radius
w0=1.31*tip_radius*tip_height/(tip_height+tip_radius+tip_radius)
wi=tip_radius/2

#refractive index of substrate (reference)
eps_star_sub=eps_sub_sub*(eps_sub_sub+(eps_sub_sub*math.tanh(q_factor*thickness_sub_sub)))/(eps_sub_sub+(eps_sub_sub*math.tanh(q_factor*thickness_sub_sub)));
rp_sub1=(eps_star_sub-1)/(eps_star_sub+1);
beta_sub=(eps_star_sub-1)/(eps_star_sub+1);

# define a bumch of zeros that can be used later

import ref
beta_sam=ref.multi(number_of_steps,q_factor,eps_top,step_define,eps_sub,eps_sub_sub,thickness_top,thickness_sub)
rp=ref.multi2(number_of_steps,q_factor,eps_top,step_define,eps_sub,eps_sub_sub,thickness_top,thickness_sub)

# Calculation of near-field
int_time=(quad(lambda x: exp(1j*x*probe_fre*harmonic), [0, time_tap]))
import near 
amp=near.part1(g,tip_distance,tip_radius,tip_height,w0,wi,C,time_tap,int_time,beta_sam,beta_sub,rp,rp_sub1,number_of_steps)
phase=near.part2(g,tip_distance,tip_radius,tip_height,w0,wi,C,time_tap,int_time,beta_sam,beta_sub,rp,rp_sub1,number_of_steps)

# plots

phase1=phase.real
for kk in range (number_of_steps):
    if phase1[kk]<=-1:
        phase1[kk]=phase1[kk]+3.14+3.14

plt.figure(1)
plt.plot (step_define,phase1)
plt.xlabel(r'$Wavenumber$ $(cm^{-1})$')
plt.ylabel(r'$Phase$ $(rad)$')
plt.title(r'$Phase$')
plt.figure(2)
plt.plot (step_define,amp)
plt.xlabel(r'$Wavenumber$ $(cm^{-1})$')
plt.ylabel(r'$Amplitude$ $(arb.unit)$')
plt.title(r'$Amplitude$')
plt.show()