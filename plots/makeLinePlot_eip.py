#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Marya Poterek
15 September 2022
Created in Python 3.7.4
This script creates contour plots for the seasonal model for chosen temperature scenarios.
"""
#%%
import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.integrate as spi
from diff_eqs_EIP import *
#%%

#%% 
def coinf_model(X, Y, eipm1, eipm2, eipm12):
    virus2introduction = max(Y,X)   
    '''overall dynamics with and without cotransmission - Fixed b'''
    #one inital mosquito with virus X (if population a million)
    INPUT = np.zeros(27)
    INPUT[1] = 1e-5
    INPUT[11] = 5e-3
    INPUT[0] = 1-sum(INPUT[1:10])
    INPUT[10] = 1-sum(INPUT[11:14])
    
    #Times
    t_start = min(X,Y); t_end = ND+min(X,Y); t_inc = TS
    #this is the time that just virus X is circulating:
    t_range1 = np.arange(t_start, virus2introduction, t_inc) 
    
    #this is the time that both viruses are circulating:
    t_range2 = np.arange(virus2introduction, t_end+t_inc, t_inc)

    params = (tempMean,tempAmp,gamma,p1H,p2H,p12H,p1M,p2M,p12M,alpha1H,alpha2H,eiph,r,eipm1, eipm2, eipm12)
    if Y != X:
    ##run the full model 
    #50% co-transmission

        RES1 = spi.odeint(diff_eqs,INPUT,t_range1,args=params)
    #Import a second infection
        INPUT2 = RES1[-1]
        INPUT2[2] = 1e-5 
        INPUT2[12] = 5e-3
        RES2 = spi.odeint(diff_eqs, INPUT2, t_range2, args = params)
        RES = np.concatenate((RES1, RES2))
        [SH, I1H, I2H, I12CH, I12SH, I1_2H, I2_1H, R1H, R2H, R12H,\
         SM, I1M, I2M, I12M, cumulative_I12CH, cumulative_I12SH, E1M, E2M,\
         E12M, E1H, E2H, E12CH, E12SH, E1_2H, E2_1H, I1E2M, I2E1M] = RES.T
    else:
        t_range = np.arange(t_start, t_end, t_inc)
        INPUT[2] = 1e-5
        INPUT[12] = 5e-3
        RES = spi.odeint(diff_eqs,INPUT,t_range, args=params)
        [SH, I1H, I2H, I12CH, I12SH, I1_2H, I2_1H, R1H, R2H, R12H,\
         SM, I1M, I2M, I12M, cumulative_I12CH, cumulative_I12SH, E1M, E2M,\
         E12M, E1H, E2H, E12CH, E12SH, E1_2H, E2_1H, I1E2M, I2E1M] = RES.T
    
    coinf_tot = (cumulative_I12SH[-1] + cumulative_I12CH[-1])
    return coinf_tot

#%%
#Define global parameters
alpha = 0.0
#These values are for Rio
tempMean = 25.1
tempAmp = 3.4

alpha1H = alpha2H = alpha1M = alpha2M = alpha
r = 5.0 
ND = 1.0*365.25
TS = 1.0
#paameters for human-> mosquito cotransmission - 60%
p1M = p2M = 0.2
p12M = 1-p1M-p2M
#parameters for mosquito -> human cotransmission - 50% initially, varied later.
p1H = p2H = 0.25
p12H = 1-p1H-p2H
eiph = 7

#%%
#To calculate gamma
g = 1./(max(quadratic(25.1,9.16,37.73,-1.48e-01),0.01))
m = (100000+35000)/(26000+83000)  #from Massad modeling paper
gamma = m*g

#%%
#To run the model and produce the contour plot

eipm_c = 2 #chikungunya
eipm_z = 14 #zika


sigma = np.linspace(0,1,11)
eipm12_range = np.zeros(len(sigma))

#Simultaneous arrival scenario
eipm1 = eipm_c #order does not matter here
eipm2 = eipm_z
X = 0
Y = 0
coinfection_out1 = np.zeros(len(sigma))

for i in range(len(sigma)):
    eipm12 = (1-sigma[i])*eipm_c + sigma[i]*eipm_z
    eipm12_range[i] = eipm12
    coinfection_out1[i] = coinf_model(X,Y,eipm1,eipm2,eipm12)


#Zika first scenario
eipm1 = eipm_z #order does matter here
eipm2 = eipm_c
X = 0
Y = 50 #50 days apart
coinfection_out2 = np.zeros(len(sigma))

for i in range(len(sigma)):
    eipm12 = (1-sigma[i])*eipm_c + sigma[i]*eipm_z
    eipm12_range[i] = eipm12
    coinfection_out2[i] = coinf_model(X,Y,eipm1,eipm2,eipm12)

#Chikungunya first scenario
eipm1 = eipm_c #order does matter here
eipm2 = eipm_z
X = 0
Y = 50 #50 days apart
coinfection_out3 = np.zeros(len(sigma))

for i in range(len(sigma)):
    eipm12 = (1-sigma[i])*eipm_c + sigma[i]*eipm_z
    eipm12_range[i] = eipm12
    coinfection_out3[i] = coinf_model(X,Y,eipm1,eipm2,eipm12)

#%%
plt.close()
plt.figure(figsize=(7,5))#,figsize=(60,40),dpi=300,pointsize=12)
plt.plot(eipm12_range, coinfection_out1, color = "red", label = "Simultaneous arrival")
plt.plot(eipm12_range, coinfection_out2, color = "blue", label = "Zika arrives first")
plt.plot(eipm12_range, coinfection_out3, color = "orange", label = "Chikungunya arrives first")
plt.axhline(y=0.083, xmin=0, xmax = 1, color='black', linestyle='dashed', linewidth=5, label ="Dengue-like parameters")
plt.legend(loc = (.55,.15))
plt.xlabel("EIP")
plt.ylabel("Cumulative incidence of co-infection")
#plt.savefig("lineplot_eip_1.png", dpi=300)
