# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 19:19:48 2023

@author: alire
"""
import math
import cmath   #import cmath for complex number operations
import numpy as np
def part1 (g,tip_distance,tip_radius,tip_height,w0,wi,C,time_tap,int_time,beta_sam,beta_sub,rp,rp_sub1,number_of_steps):

    f01=(g-((2*tip_distance+w0+tip_radius)/2/tip_height))
    f02=math.log(4*tip_height/(4*tip_distance+2*w0+tip_radius))
    f03=math.log(4*tip_height/tip_radius)
    f0=f01*f02/f03
    fh1=(g-((2*tip_distance+wi+tip_radius)/2/tip_height))
    fh2=math.log(4*tip_height/(4*tip_distance+2*tip_radius))
    fh3=math.log(4*tip_height/tip_radius)
    fh=fh1*fh2/fh3
    xx_new_zero=np.zeros(number_of_steps)
    alfaefective_sam=xx_new_zero+1j*xx_new_zero; alfaefective=xx_new_zero+1j*xx_new_zero
    alfaefective_sub1=xx_new_zero+1j*xx_new_zero; alfaefective_sub=xx_new_zero+1j*xx_new_zero
    top_alfa=xx_new_zero+1j*xx_new_zero; sub_alfa=xx_new_zero+1j*xx_new_zero
    phase=xx_new_zero+1j*xx_new_zero; amp=xx_new_zero+1j*xx_new_zero
    for kk in range (number_of_steps):
        alfaefective_sam[kk]=C*(1+(f0*beta_sam[kk]/2/(1-fh*beta_sam[kk])));
        alfaefective[kk]=1/time_tap*int_time*alfaefective_sam[kk];
        alfaefective_sub1[kk]=C*(1+(f0*beta_sub/2/(1-fh*beta_sub))); 
        alfaefective_sub[kk]=1/time_tap*int_time*alfaefective_sub1[kk];
        top_alfa[kk]=alfaefective[kk]*(1+rp[kk])*(1+rp[kk])
        sub_alfa[kk]=alfaefective_sub[kk]*(1+rp_sub1)*(1+rp_sub1)
        phase[kk]=cmath.phase(top_alfa[kk]/sub_alfa[kk])
        amp[kk]=abs(top_alfa[kk]/sub_alfa[kk])
    return (amp)
def part2 (g,tip_distance,tip_radius,tip_height,w0,wi,C,time_tap,int_time,beta_sam,beta_sub,rp,rp_sub1,number_of_steps):

    f01=(g-((2*tip_distance+w0+tip_radius)/2/tip_height))
    f02=math.log(4*tip_height/(4*tip_distance+2*w0+tip_radius))
    f03=math.log(4*tip_height/tip_radius)
    f0=f01*f02/f03
    fh1=(g-((2*tip_distance+wi+tip_radius)/2/tip_height))
    fh2=math.log(4*tip_height/(4*tip_distance+2*tip_radius))
    fh3=math.log(4*tip_height/tip_radius)
    fh=fh1*fh2/fh3
    xx_new_zero=np.zeros(number_of_steps)
    alfaefective_sam=xx_new_zero+1j*xx_new_zero; alfaefective=xx_new_zero+1j*xx_new_zero
    alfaefective_sub1=xx_new_zero+1j*xx_new_zero; alfaefective_sub=xx_new_zero+1j*xx_new_zero
    top_alfa=xx_new_zero+1j*xx_new_zero; sub_alfa=xx_new_zero+1j*xx_new_zero
    phase=xx_new_zero+1j*xx_new_zero; amp=xx_new_zero+1j*xx_new_zero
    for kk in range (number_of_steps):
        alfaefective_sam[kk]=C*(1+(f0*beta_sam[kk]/2/(1-fh*beta_sam[kk])));
        alfaefective[kk]=1/time_tap*int_time*alfaefective_sam[kk];
        alfaefective_sub1[kk]=C*(1+(f0*beta_sub/2/(1-fh*beta_sub))); 
        alfaefective_sub[kk]=1/time_tap*int_time*alfaefective_sub1[kk];
        top_alfa[kk]=alfaefective[kk]*(1+rp[kk])*(1+rp[kk])
        sub_alfa[kk]=alfaefective_sub[kk]*(1+rp_sub1)*(1+rp_sub1)
        phase[kk]=cmath.phase(top_alfa[kk]/sub_alfa[kk])
        amp[kk]=abs(top_alfa[kk]/sub_alfa[kk])
    return (phase)
