#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: mlpoterek
"""
import scipy.integrate as spi
import numpy as np
import sys

def coinf_model(param_values):
    Y = np.zeros((param_values.shape[0], 5))
    for i, X in enumerate(param_values):
        m = X[0]
        a = X[1]
        b1 = b2 = b12 = X[2]
        c1 = c2 = c12 = X[3]
        r = X[4]
        g = X[5]
        eiph = X[6]
        eipm = X[7]
        virus2introduction = X[8]
        alpha1H = alpha2H = alpha1M = alpha2M = X[9]
        ND = 1.75*365.25
        TS = 1.0
        
        def diff_eqs(INP,t):  
            '''The main set of equations'''
            [SH, I1H, I2H, I12CH, I12SH, I1_2H, I2_1H, R1H, R2H, R12H,\
             SM, I1M, I2M, I12M, cumulative_I12CH, cumulative_I12SH, E1M, E2M,\
             E12M, E1H, E2H, E12CH, E12SH, E1_2H, E2_1H, I1E2M, I2E1M] = INP
             
            lambda1H = m*a*b1*(I1M+I12M*p1H)
            lambda2H = m*a*b2*(I2M+I12M*p2H)
            lambda12H = m*a*b12*I12M*p12H
             
            lambda1M = a*c1*(I1H+I2_1H+(I12CH+I12SH)*p1M)
            lambda2M = a*c2*(I2H+I1_2H+(I12CH+I12SH)*p2M)
            lambda12M = a*c12*(I12CH+I12SH)*p12M
             
            dSH = -(lambda1H+lambda2H+lambda12H)*SH
            dE1H = lambda1H*SH - ((lambda2H + lambda12H)*(1-alpha1H) + (1./eiph))*E1H  
            dE2H = lambda2H*SH - ((lambda1H + lambda12H)*(1-alpha2H) + (1./eiph))*E2H  
            dE12CH = lambda12H*SH + lambda12H*((1-alpha2H)*I2H + (1-alpha1H)*I1H) - (1./eiph)*E12CH
            dE12SH = lambda2H*(1-alpha1H)*I1H + lambda1H*(1-alpha2H)*I2H - (1./eiph)*E12SH
            dI1H = (1./eiph)*E1H - (1./r)*I1H
            dI2H = (1./eiph)*E2H - (1./r)*I2H
            dI12CH = (1./eiph)*E12CH - (1./r)*I12CH
            dI12SH = (1./eiph)*E12SH - (1./r)*I12SH
            dR1H = (1./r)*I1H - (lambda12H+lambda2H)*(1-alpha1H)*R1H
            dR2H = (1./r)*I2H - (lambda12H+lambda1H)*(1-alpha2H)*R2H
            dR12H = (1./r)*(I2_1H + I1_2H + I12CH + I12SH)
            dE1_2H = (lambda12H+lambda2H)*(1-alpha1H)*R1H - (1./eiph)*E1_2H
            dE2_1H = (lambda12H+lambda1H)*(1-alpha2H)*R2H - (1./eiph)*E2_1H
            dI1_2H = (1./eiph)*E1_2H - (1./r)*I1_2H
            dI2_1H = (1./eiph)*E2_1H - (1./r)*I2_1H
             
             #Mosquito
            dSM = g*(I1M+I2M+I12M) - (lambda1M+lambda2M+lambda12M)*SM
            dE1M = lambda1M*SM - ((lambda2M+lambda12M)*(1-alpha1M) + g)*E1M - (1./eipm)*E1M
            dE2M = lambda2M*SM - ((lambda1M+lambda12M)*(1-alpha2M) + g)*E2M - (1./eipm)*E2M
            dE12M = lambda12M*SM + (lambda2M+lambda12M)*(1-alpha1M)*E1M + (lambda1M+lambda12M)*(1-alpha2M)*E2M - g*E12M - (1./eipm)*E12M 
            dI1M = (1./eipm)*E1M - lambda2M*I1M - g*I1M 
            dI2M = (1./eipm)*E2M - lambda1M*I2M - g*I2M
            dI12M = (1./eipm)*(E12M + I1E2M + I2E1M) - g*I12M
            dI1E2M = lambda2M*I1M - (1./eipm)*I1E2M - g*I1E2M
            dI2E1M = lambda1M*I2M - (1./eipm)*I2E1M - g*I2E1M
            
            dcumulative_I12CH = lambda12H*SH + lambda12H*((1-alpha2H)*I2H + (1-alpha1H)*I1H)
            dcumulative_I12SH = lambda2H*(1-alpha1H)*I1H + lambda1H*(1-alpha2H)*I2H
    
            Y = np.array([dSH, dI1H, dI2H, dI12CH, dI12SH, dI1_2H, dI2_1H, dR1H, dR2H,\
                       dR12H, dSM, dI1M, dI2M, dI12M, dcumulative_I12CH, dcumulative_I12SH, dE1M, dE2M, \
                       dE12M, dE1H, dE2H, dE12CH, dE12SH, dE1_2H, dE2_1H, dI1E2M, dI2E1M])
    
            return Y

        '''overall dynamics with and without cotransmission - Fixed b'''
        #one inital mosquito with virus X (if population a million)
        INPUT = np.zeros(27)
        INPUT[1] = 1e-6
        INPUT[0] = 1-sum(INPUT[1:10])
        INPUT[10] = 1-sum(INPUT[11:14])

    #parameters for human-> mosquito cotransmission - 60%
	p1M = p2M = 0.2
        p12M = 1-p1M-p2M
        #parameters for mosquito -> human cotransmission - 50% initially, varied later.
        p1H = p2H = 0.25
        p12H = 1-p1H-p2H
        
        t_start = 0.0; t_end = ND; t_inc = TS
        t_range1 = np.arange(t_start, virus2introduction, t_inc) 
        t_range2 = np.arange(virus2introduction, t_end+t_inc, t_inc)
        t_range = np.concatenate((t_range1, t_range2))
                
        #50% co-transmission
        RES1 = spi.odeint(diff_eqs,INPUT,t_range)
        [SH, I1H, I2H, I12CH, I12SH, I1_2H, I2_1H, R1H, R2H, R12H,\
         SM, I1M, I2M, I12M, cumulative_I12CH, cumulative_I12SH, E1M, E2M,\
         E12M, E1H, E2H, E12CH, E12SH, E1_2H, E2_1H, I1E2M, I2E1M] = RES1.T
                 
        #50% co-transmission
        RES1 = spi.odeint(diff_eqs,INPUT,t_range1)    
        #Import a second infection
        INPUT2 = RES1[-1] 
        INPUT2[2] = 1e-6 
        RES2 = spi.odeint(diff_eqs, INPUT2, t_range2)
        RES = np.concatenate((RES1, RES2))
        [SH, I1H, I2H, I12CH, I12SH, I1_2H, I2_1H, R1H, R2H, R12H,\
         SM, I1M, I2M, I12M, cumulative_I12CH, cumulative_I12SH, E1M, E2M,\
         E12M, E1H, E2H, E12CH, E12SH, E1_2H, E2_1H, I1E2M, I2E1M] = RES.T
         
        #total number infected w 1 
        inf_1 = R1H[-1]*(1e6)
        #total number infected w 2 
        inf_2 = R2H[-1]*(1e6)
        #total sequential coinfections 
        seq_inf = cumulative_I12SH[-1]*(1e6)
        #total cotransmitted coinfections 
        cot_inf = cumulative_I12CH[-1]*(1e6)
        #proportion of all infections that are coinfectiodns
        coinf_prop = (cumulative_I12SH[-1] + cumulative_I12CH[-1])/(R12H[-1] + R1H[-1] + R2H[-1])
        output_vec = [inf_1, inf_2, seq_inf, cot_inf, coinf_prop]
        Y[i] = np.array(output_vec) 
    	return Y
param_values = np.loadtxt("file"+sys.argv[1]+".txt", float)
Z = coinf_model(param_values)
np.savetxt(sys.argv[1]+".out", Z)
