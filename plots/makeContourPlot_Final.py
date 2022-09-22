#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Marya Poterek
21 February 2021
Created in Python 3.7.4
This script creates contour plots for the seasonal model for chosen temperature scenarios.
"""
#%%
import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.integrate as spi
from diff_eqs import *
#%%

#%% 
def coinf_model(X, Y):
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

    params = (tempMean,tempAmp,gamma,p1H,p2H,p12H,p1M,p2M,p12M,alpha1H,alpha2H,eiph,r)
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
    else: #I think this gives the same amount of co-infection
        #print("Y",Y,"X",X)
        t_range = np.arange(t_start, t_end, t_inc)
        INPUT[2] = 1e-5
        INPUT[12] = 5e-3
        RES = spi.odeint(diff_eqs,INPUT,t_range, args=params)
        [SH, I1H, I2H, I12CH, I12SH, I1_2H, I2_1H, R1H, R2H, R12H,\
         SM, I1M, I2M, I12M, cumulative_I12CH, cumulative_I12SH, E1M, E2M,\
         E12M, E1H, E2H, E12CH, E12SH, E1_2H, E2_1H, I1E2M, I2E1M] = RES.T
    
    coinf_tot = (cumulative_I12SH[-1] + cumulative_I12CH[-1])
    return coinf_tot

#Define global parameters
alpha = 0.0

#Values for Rio de Janeiro
#tempMean = 25.1
#tempAmp = 3.4

#Values are for Beijing
tempMean = 13.8
tempAmp = 14.7
alpha1H = alpha2H = alpha1M = alpha2M = alpha
r = 5.0 
ND = 1.0*365.25
TS = 1.0
#parameters for human-> mosquito cotransmission - 60%
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
length = 365
virusAintro = np.zeros(length)
virusBintro = np.zeros(length)
for day in range(length):
    virusAintro[day] = day
    virusBintro[day] = day
virusAintro_grid, virusBintro_grid = np.meshgrid(virusAintro,virusBintro)
coinfection_out = np.zeros((length*length)).reshape(length,length)
for i in range(len(virusAintro_grid)):
    for j in range(len(virusBintro_grid)):
        if j >= i:
            print("here")
            coinfection_out[i][j] = coinf_model(virusBintro_grid[i][j],virusAintro_grid[i][j])
            coinfection_out[j][i] = coinfection_out[i][j]

fig, axes = plt.subplots(nrows=2, ncols=3,
                         gridspec_kw={'height_ratios': [7,1], 
                                      'width_ratios': [1,7,1]},
                         figsize=(0.7*9,0.7*8))
#Lines corresponding to slice plot
x = np.linspace(0, 362, 365)
y = np.linspace(0, 362, 365)
x1 = np.linspace(0, 337,365)
y1 = np.linspace(25,362,365)
x2 = np.linspace(0, 312, 365)
y2 = np.linspace(50, 362, 365)
contourplot = plt.contourf(virusAintro_grid,virusBintro_grid,coinfection_out,40,cmap = "viridis")
plt.gca().set_visible(False)
axes[0,1].contourf(virusAintro_grid,virusBintro_grid,coinfection_out,40,cmap = "viridis")
axes[0,1].plot(x, y, c = "red")
axes[0,1].plot(x1,y1, c = "orange")
axes[0,1].plot(x2,y2, c = "green")
fig.colorbar(contourplot,cax=axes[0][2]).set_label('Total annual co-infections', rotation=270, labelpad = 15)
axes[1,0].axis('off')
axes[0,1].axis('off')
axes[1,2].axis('off')
tvec = np.arange(1,366)
axes[1,1].plot(tvec, temp_func(tempMean,tempAmp,tvec))
axes[1,1].set_xlim([0,365])
axes[1,1].spines['top'].set_visible(False)
axes[1,1].spines['right'].set_visible(False)
axes[1,1].spines['bottom'].set_visible(False)
axes[1,1].spines['left'].set_visible(False)
axes[1,1].set_ylabel("T (C)")
axes[1,1].set_xlabel("Day of arrival")
axes[0,0].plot(temp_func(tempMean,tempAmp,tvec), tvec)
axes[0,0].set_ylim([0,365])
axes[0,0].spines['top'].set_visible(False)
axes[0,0].spines['right'].set_visible(False)
axes[0,0].spines['bottom'].set_visible(False)
axes[0,0].spines['left'].set_visible(False)
axes[0,0].set_xlabel("T (C)")
axes[0,0].set_ylabel("Day of arrival")

#plt.savefig("figure1_A.png", dpi=300) #rio
#plt.savefig("figure1_B.png",dpi=300) #beijing
