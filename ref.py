# -*- coding: utf-8 -*-
"""
Created on Sun Jan 22 17:49:58 2023

@author: alire
"""
"""
This functin use the dielectric function of 3 surfaces and thicknesses to calculate
the near-field and far-field reflectivity of surface
"""
import math
import cmath   #import cmath for complex number operations
import numpy as np
"""
Summery:this function produce the near field (beta) and farfield (rp) refractive of the sample


"""
def multi (number_of_steps,q_factor,eps_top,step_define,eps_sub,eps_sub_sub,thickness_top,thickness_sub):
    xx_new_zero=np.zeros(number_of_steps)
    b12 = xx_new_zero+1j*xx_new_zero; b23 = xx_new_zero+1j*xx_new_zero; b34 = xx_new_zero+1j*xx_new_zero
    ba = xx_new_zero+1j*xx_new_zero; bb = xx_new_zero+1j*xx_new_zero
    beta_sam = xx_new_zero+1j*xx_new_zero

    for counting in range (number_of_steps):
        b12[counting]=(eps_top[counting]-1)/(eps_top[counting]+1)
        b23[counting]=(eps_sub[counting]-eps_top[counting])/(eps_sub[counting]+eps_top[counting])
        b34[counting]=(eps_sub_sub-eps_sub[counting])/(eps_sub_sub+eps_sub[counting])
        ba[counting]=(b12[counting]*b23[counting]*b34[counting]*math.exp(-2*q_factor*(thickness_sub)))+(b23[counting]*math.exp(-2*q_factor*(thickness_top)))+(b34[counting]*math.exp(-2*q_factor*(thickness_top+thickness_sub)))
        bb[counting]=(b12[counting]*b23[counting]*math.exp(-2*q_factor*(thickness_top)))+(b12[counting]*b34[counting]*math.exp(-2*q_factor*(thickness_top+thickness_sub)))+(b23[counting]*b34[counting]*math.exp(-2*q_factor*(thickness_sub)))
        beta_sam[counting]=(b12[counting]+ba[counting])/(1+bb[counting])
    return (beta_sam)
# refractive index of top material

def multi2 (number_of_steps,q_factor,eps_top,step_define,eps_sub,eps_sub_sub,thickness_top,thickness_sub):
    xx_new_zero=np.zeros(number_of_steps)
    k1 = xx_new_zero+1j*xx_new_zero; k2 = xx_new_zero+1j*xx_new_zero; k3 = xx_new_zero+1j*xx_new_zero; k4 = xx_new_zero+1j*xx_new_zero
    r12 = xx_new_zero+1j*xx_new_zero; r23 = xx_new_zero+1j*xx_new_zero; r34 = xx_new_zero+1j*xx_new_zero
    r0a = xx_new_zero+1j*xx_new_zero; r1a = xx_new_zero+1j*xx_new_zero; r2a = xx_new_zero+1j*xx_new_zero
    ra = xx_new_zero+1j*xx_new_zero; rb = xx_new_zero+1j*xx_new_zero
    rb1 = xx_new_zero+1j*xx_new_zero; rb2 = xx_new_zero+1j*xx_new_zero; rb3 = xx_new_zero+1j*xx_new_zero
    rp = xx_new_zero+1j*xx_new_zero; 
    for counting in range (number_of_steps):
        k1[counting]=cmath.sqrt(step_define[counting]*step_define[counting]/3e8/3e8-q_factor*q_factor)
        k2[counting]=cmath.sqrt(eps_top[counting]*step_define[counting]*step_define[counting]/3e8/3e8-q_factor*q_factor)
        k3[counting]=cmath.sqrt(eps_sub[counting]*step_define[counting]*step_define[counting]/3e8/3e8-q_factor*q_factor)
        k4[counting]=cmath.sqrt(eps_sub_sub*step_define[counting]*step_define[counting]/3e8/3e8-q_factor*q_factor)
        r12[counting]=(eps_top[counting]*k1[counting]-k2[counting])/(eps_top[counting]*k1[counting]+k2[counting])
        r23[counting]=(eps_sub[counting]*k2[counting]-eps_top[counting]*k3[counting])/(eps_sub[counting]*k2[counting]+eps_top[counting]*k3[counting])
        r34[counting]=(eps_sub_sub*k3[counting]-eps_sub[counting]*k4[counting])/(eps_sub_sub*k3[counting]+eps_sub[counting]*k4[counting])
        r0a[counting]=cmath.exp(k2[counting]*1j*2*thickness_sub)
        r1a[counting]=r23[counting]*cmath.exp(1j*2*k2[counting]*thickness_top)
        r2a[counting]=r34[counting]*cmath.exp(1j*2*k2[counting]*(thickness_top+thickness_sub))
        ra[counting]=r12[counting]*r23[counting]*r34[counting]*r0a[counting]+r1a[counting]+r2a[counting]
        rb1[counting]=r12[counting]*r23[counting]*cmath.exp(1j*2*k2[counting]*thickness_top)
        rb2[counting]=r12[counting]*r34[counting]*cmath.exp(1j*2*k2[counting]*(thickness_top+thickness_sub))
        rb3[counting]=r23[counting]*r34[counting]*cmath.exp(1j*2*k2[counting]*(thickness_sub))
        rb[counting]=rb1[counting]+rb2[counting]+rb3[counting]
        rp[counting]=(r12[counting]+ra[counting])/(1+rb[counting])
    return (rp)