#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: mlpoterek
"""
#%%
import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.integrate as spi
from diff_eqs import *

#%% 
def coinf_model(start):
    INPUT = np.zeros(27)
    INPUT[1] = 1e-5
    INPUT[11] = 5e-3
    INPUT[0] = 1-sum(INPUT[1:10])
    INPUT[10] = 1-sum(INPUT[11:14])

    #Times
    t_start = start; t_end = ND; t_inc = TS
    #this is the time that just virus X is circulating:
    t_range1 = np.arange(t_start, t_start+virus2introduction, t_inc) 
#this is the time that both viruses are circulating:
    t_range2 = np.arange(t_start+virus2introduction, t_end+t_inc, t_inc)
    #t_range = np.concatenate((t_range1, t_range2))
    
#50% co-transmission
    params = (tempMean,tempAmp,gamma,p1H,p2H,p12H,p1M,p2M,p12M,alpha1H,alpha2H,eiph,r)
    RES1 = spi.odeint(diff_eqs,INPUT,t_range1,args=params)

#Import a second infection
    INPUT2 = RES1[-1] 
    INPUT2[2] = 1e-5 
    INPUT2[12] = 5e-3
    RES2 = spi.odeint(diff_eqs, INPUT2, t_range2,args=params)
    RES = np.concatenate((RES1, RES2))
    [SH, I1H, I2H, I12CH, I12SH, I1_2H, I2_1H, R1H, R2H, R12H,\
     SM, I1M, I2M, I12M, cumulative_I12CH, cumulative_I12SH, E1M, E2M,\
     E12M, E1H, E2H, E12CH, E12SH, E1_2H, E2_1H, I1E2M, I2E1M] = RES.T
    
    
    #proportion of all infections that are coinfections
    inf1 = R1H[-1]+I1_2H[-1]+R12H[-1]
    inf2 = R2H[-1]+I2_1H[-1]+R12H[-1]
    seq_tot = cumulative_I12SH[-1]
    cot_tot = cumulative_I12CH[-1]
    coinf_tot = (cumulative_I12SH[-1] + cumulative_I12CH[-1])
    coinf_prop = (cumulative_I12SH[-1] + cumulative_I12CH[-1])/(R12H[-1] + R1H[-1] + R2H[-1])
    cot_seq_ratio = cot_tot / seq_tot
    return inf1, coinf_prop, cot_seq_ratio
    
#%%
#Define temp-dependent global parameters
alpha = 0.0
tempMean = 25.1
tempAmp = 0
alpha1H = alpha2H = alpha
r = 5.0 
ND = 1.0*365.25
TS = 1.0
#parameters for human-> mosquito cotransmission - 60%
p1M = p2M = 0.2
p12M = 1-p1M-p2M
#parameters for mosquito -> human cotransmission - 50% initially, varied later.
eiph = 7
virus2introduction = 30

#To calculate gamma
g = 1./(max(quadratic(25.1,9.16,37.73,-1.48e-01),0.01))
m = (100000+35000)/(26000+83000)  #from Massad modeling paper
gamma = m*g

#To calculate the total co-infection observed using the baseline parameters
#Non seasonal model with start date = 0


#Reset values after generating gamma to explore individual parameter spaces
tempAmp = 3.4
p1H = p2H = 0.25
p12H = 1-p1H-p2H


#%%
#To create a 3x3 figure
#alpha
o1b=np.zeros(len(np.linspace(0,1,21)))
o2b=np.zeros(len(np.linspace(0,1,21)))
o3b=np.zeros(len(np.linspace(0,1,21)))
count = 0
for a in np.linspace(0,1,21):
    alpha = a
    alpha1H = alpha2H = alpha1M = alpha2M = alpha
    o1b[count],o2b[count],o3b[count] = coinf_model(0.0)
    count = count+1

    
o1c=np.zeros(len(np.linspace(0,1,21)))
o2c=np.zeros(len(np.linspace(0,1,21)))
o3c=np.zeros(len(np.linspace(0,1,21)))
count = 0
for a in np.linspace(0,1,21):
    alpha = a
    alpha1H = alpha2H = alpha1M = alpha2M = alpha
    o1c[count],o2c[count],o3c[count] = coinf_model(ND*.25)
    count = count+1

    
o1d=np.zeros(len(np.linspace(0,1,21)))
o2d=np.zeros(len(np.linspace(0,1,21)))
o3d=np.zeros(len(np.linspace(0,1,21)))
count = 0
for a in np.linspace(0,1,21):
    alpha = a
    alpha1H = alpha2H = alpha1M = alpha2M = alpha
    o1d[count],o2d[count],o3d[count] = coinf_model(ND*.5)
    count = count+1

o1e=np.zeros(len(np.linspace(0,1,21)))
o2e=np.zeros(len(np.linspace(0,1,21)))
o3e=np.zeros(len(np.linspace(0,1,21)))
count = 0
for a in np.linspace(0,1,21):
    alpha = a
    alpha1H = alpha2H = alpha1M = alpha2M = alpha
    o1e[count],o2e[count],o3e[count] = coinf_model(.75*ND)
    count = count+1
 
 #r
alpha = alpha1H = alpha2H = alpha1M = alpha2M = alpha= 0.0
o4b=np.zeros(len(np.linspace(3,10,8)))
o5b=np.zeros(len(np.linspace(3,10,8)))
o6b=np.zeros(len(np.linspace(3,10,8)))
count = 0
for a in np.linspace(3,10,8):
    r = a
    o4b[count],o5b[count],o6b[count] = coinf_model(0.0)
    count = count+1

o4c=np.zeros(len(np.linspace(3,10,8)))
o5c=np.zeros(len(np.linspace(3,10,8)))
o6c=np.zeros(len(np.linspace(3,10,8)))
count = 0
for a in np.linspace(3,10,8):
    r = a
    o4c[count],o5c[count],o6c[count] = coinf_model(.25*ND)
    count = count+1

o4d=np.zeros(len(np.linspace(3,10,8)))
o5d=np.zeros(len(np.linspace(3,10,8)))
o6d=np.zeros(len(np.linspace(3,10,8)))
count = 0
for a in np.linspace(3,10,8):
    r = a
    o4d[count],o5d[count],o6d[count] = coinf_model(.5*ND)
    count = count+1
    
o4e=np.zeros(len(np.linspace(3,10,8)))
o5e=np.zeros(len(np.linspace(3,10,8)))
o6e=np.zeros(len(np.linspace(3,10,8)))
count = 0
for a in np.linspace(3,10,8):
    r = a
    o4e[count],o5e[count],o6e[count] = coinf_model(.75*ND)
    count = count+1
    
#eiph
r = 5
o7b=np.zeros(len(np.linspace(4,8,5)))
o8b=np.zeros(len(np.linspace(4,8,5)))
o9b=np.zeros(len(np.linspace(4,8,5)))
count = 0
for a in np.linspace(4,8,5):
    eiph = a
    o7b[count],o8b[count],o9b[count] = coinf_model(0.0)
    count = count+1

o7c=np.zeros(len(np.linspace(4,8,5)))
o8c=np.zeros(len(np.linspace(4,8,5)))
o9c=np.zeros(len(np.linspace(4,8,5)))
count = 0
for a in np.linspace(4,8,5):
    eiph = a
    o7c[count],o8c[count],o9c[count] = coinf_model(.25*ND)
    count = count+1

o7d=np.zeros(len(np.linspace(4,8,5)))
o8d=np.zeros(len(np.linspace(4,8,5)))
o9d=np.zeros(len(np.linspace(4,8,5))) 
count = 0
for a in np.linspace(4,8,5):
    eiph = a
    o7d[count],o8d[count],o9d[count] = coinf_model(.5*ND)
    count = count+1

o7e=np.zeros(len(np.linspace(4,8,5)))
o8e=np.zeros(len(np.linspace(4,8,5)))
o9e=np.zeros(len(np.linspace(4,8,5)))
count = 0
for a in np.linspace(4,8,5):
    eiph = a
    o7e[count],o8e[count],o9e[count] = coinf_model(.75*ND)
    count = count+1
    
#%%
#Define temp-independent global parameters
alpha = 0.0
tempMean = 25.1
tempAmp = 0.0
alpha1H = alpha2H= alpha
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
virus2introduction = 30

#To calculate gamma
g = 1./(max(quadratic(24.3,9.16,37.73,-1.48e-01),0.01))
m = (100000+35000)/(26000+83000)  #from Massad modeling paper
gamma = m*g
#%%

#alpha
o1=np.zeros(len(np.linspace(0,1,21)))
o2=np.zeros(len(np.linspace(0,1,21)))
o3=np.zeros(len(np.linspace(0,1,21)))

count = 0
for a in np.linspace(0,1,21):
    alpha = a
    alpha1H = alpha2H = alpha1M = alpha2M = alpha
    o1[count],o2[count],o3[count] = coinf_model(0.0)
    count = count+1
 
 #r
alpha = alpha1H = alpha2H = alpha1M = alpha2M = alpha= 0.0
o4=np.zeros(len(np.linspace(3,10,8)))
o5=np.zeros(len(np.linspace(3,10,8)))
o6=np.zeros(len(np.linspace(3,10,8)))

count = 0
for a in np.linspace(3,10,8):
    r = a
    o4[count],o5[count],o6[count] = coinf_model(0.0)
    count = count+1
    
#eiph
r = 5
o7=np.zeros(len(np.linspace(4,8,5)))
o8=np.zeros(len(np.linspace(4,8,5)))
o9=np.zeros(len(np.linspace(4,8,5)))
 
count = 0
for a in np.linspace(4,8,5):
    eiph = a
    o7[count],o8[count],o9[count] = coinf_model(0.0)
    count = count+1
#%% 
#to plot all     
fig, axs  = plt.subplots(3,3, figsize=(12,8))#, sharex = 'all', sharey = "none")
axs[0, 0].plot(np.linspace(0,1,21),o1, color = "black", linestyle = "dashed")
axs[0, 0].plot(np.linspace(0,1,21),o1b, color = "black")
axs[0, 0].plot(np.linspace(0,1,21),o1c, color = "dimgray")
axs[0, 0].plot(np.linspace(0,1,21),o1d, color = "darkgray")
axs[0, 0].plot(np.linspace(0,1,21),o1e, color = "lightgray")
axs[0,0].axvline(0.0, color= "gray")
axs[0,0].ticklabel_format(useOffset=False)
axs[0,0].set_ylabel('Virus A Infections', rotation=90,va ='bottom', fontsize = 12)# fontsize=20, labelpad=20)
axs[0,0].set_title('A', fontsize=12, loc="left")
axs[0,0].set_xticks([])

axs[0, 1].plot(np.linspace(3,10,8),o4, color = "black", linestyle = "dashed")
axs[0, 1].plot(np.linspace(3,10,8),o4b, color = "black")
axs[0, 1].plot(np.linspace(3,10,8),o4c, color = "dimgray")
axs[0, 1].plot(np.linspace(3,10,8),o4d, color = "darkgray")
axs[0, 1].plot(np.linspace(3,10,8),o4e, color = "lightgray")
axs[0,1].axvline(5, color= "gray")
axs[0,1].ticklabel_format(useOffset=False)
axs[0,1].set_title('B',fontsize=12, loc="left")
axs[0,1].set_xticks([])

axs[0, 2].plot(np.linspace(4,8,5),o7, color = "black", linestyle = "dashed")
axs[0, 2].plot(np.linspace(4,8,5),o7b,color = "black")
axs[0, 2].plot(np.linspace(4,8,5),o7c,color = "dimgray")
axs[0, 2].plot(np.linspace(4,8,5),o7d,color = "darkgray")
axs[0, 2].plot(np.linspace(4,8,5),o7e,color = "lightgray")
axs[0,2].axvline(7, color= "gray")
axs[0, 2].ticklabel_format(useOffset=False)
axs[0,2].set_title('C',fontsize=12,loc="left")
axs[0,2].set_xticks([])

axs[1, 0].plot(np.linspace(0,1,21),o2, color = "black", linestyle = "dashed")
axs[1, 0].plot(np.linspace(0,1,21),o2b, color = "black")
axs[1, 0].plot(np.linspace(0,1,21),o2c, color = "dimgray")
axs[1, 0].plot(np.linspace(0,1,21),o2d, color = "darkgray")
axs[1, 0].plot(np.linspace(0,1,21),o2e, color = "lightgray")
axs[1,0].axvline(0.0, color= "gray")
axs[1, 0].ticklabel_format(useOffset=False)
axs[1,0].set_ylabel('Proportion of co-infections', rotation=90,va ='bottom', fontsize = 12)# fontsize=20, labelpad=20)
axs[1,0].set_title('D',fontsize=12,loc="left")
axs[1,0].set_xticks([])

axs[1, 1].plot(np.linspace(3,10,8),o5, color = "black", linestyle = "dashed")
axs[1, 1].plot(np.linspace(3,10,8),o5b, color = "black")
axs[1, 1].plot(np.linspace(3,10,8),o5c, color = "dimgray")
axs[1, 1].plot(np.linspace(3,10,8),o5d, color = "darkgray")
axs[1, 1].plot(np.linspace(3,10,8),o5e, color = "lightgray")
axs[1,1].axvline(5, color= "gray")
axs[1, 1].ticklabel_format(useOffset=False)
axs[1,1].set_title('E',fontsize=12,loc="left")
axs[1,1].set_xticks([])

axs[1, 2].plot(np.linspace(4,8,5),o8, color = "black", linestyle = "dashed")
axs[1, 2].plot(np.linspace(4,8,5),o8b, color = "black")
axs[1, 2].plot(np.linspace(4,8,5),o8c, color = "dimgray")
axs[1, 2].plot(np.linspace(4,8,5),o8d, color = "darkgray")
axs[1, 2].plot(np.linspace(4,8,5),o8e, color = "lightgray")
axs[1,2].axvline(7, color= "gray")
axs[1, 2].ticklabel_format(useOffset=False)
axs[1,2].set_title('F',fontsize=12,loc="left")
axs[1,2].set_xticks([])

axs[2, 0].plot(np.linspace(0,1,21),o3, color = "black", linestyle = "dashed", label = "Nonseasonal model")
axs[2, 0].plot(np.linspace(0,1,21),o3b, color = "black", label = "Seasonal, start = day 0")
axs[2, 0].plot(np.linspace(0,1,21),o3c, color = "dimgray", label = "Seasonal, start = day 90")
axs[2, 0].plot(np.linspace(0,1,21),o3d, color = "darkgray", label = "Seasonal, start = day 180")
axs[2, 0].plot(np.linspace(0,1,21),o3e, color = "lightgray", label = "Seasonal, start = day 270")
axs[2,0].axvline(0.0, color= "gray")
axs[2,0].legend(loc="upper left", fontsize = "small")
axs[2,0].ticklabel_format(useOffset=False)
axs[2,0].set_ylabel('Cot:Seq co-infections', rotation=90,va ='bottom', fontsize = 12)# fontsize=20, labelpad=20)
axs[2,0].set_title('G',fontsize=12,loc="left")
axs[2, 0].set_xlabel('Cross protection (α)', fontsize =12)

axs[2, 1].plot(np.linspace(3,10,8),o6,color = "black", linestyle = "dashed")
axs[2, 1].plot(np.linspace(3,10,8),o6b,color = "black")
axs[2, 1].plot(np.linspace(3,10,8),o6c,color = "dimgray")
axs[2, 1].plot(np.linspace(3,10,8),o6d,color = "darkgray")
axs[2, 1].plot(np.linspace(3,10,8),o6e,color = "lightgray")
axs[2,1].axvline(5, color= "gray")
axs[2,1].ticklabel_format(useOffset=False)
axs[2,1].set_title('H',fontsize=12,loc="left")
axs[2, 1].set_xlabel('Infectious period length (r)', fontsize = 12)

axs[2, 2].plot(np.linspace(4,8,5),o9, color = "black", linestyle = "dashed")
axs[2, 2].plot(np.linspace(4,8,5),o9b, color = "black")
axs[2, 2].plot(np.linspace(4,8,5),o9c, color = "dimgray")
axs[2, 2].plot(np.linspace(4,8,5),o9d, color = "darkgray")
axs[2, 2].plot(np.linspace(4,8,5),o9e, color = "lightgray")
axs[2,2].axvline(7, color= "gray")
axs[2, 2].ticklabel_format(useOffset=False)
axs[2,2].set_title("J",fontsize=12,loc="left")
axs[2, 2].set_xlabel('Intrinsic incubation period (ε)', fontsize= 12)
