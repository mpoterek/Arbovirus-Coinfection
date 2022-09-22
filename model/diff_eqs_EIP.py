import math
import numpy as np
def temp_func(temp_mean, amp,t):
    T = 0.5*(amp)*np.sin(t*(2*math.pi)/365) + temp_mean
    return T
def briere(T,Tmin,Tmax,c):
    if (T<Tmin or T>Tmax):
        Bt = 0
    else:
   #     Bt = c*T*(T-Tmin)*math.sqrt(Tmax-T)
        Bt = max(c*T*(T-Tmin)*math.sqrt(Tmax-T),0) 
    return Bt
def quadratic(T,Tmin,Tmax,c):
    Q = max(c*(T-Tmax)*(T-Tmin),0)
    return Q
def diff_eqs(INP,t,tempMean,tempAmp,gamma,p1H,p2H,p12H,p1M,p2M,p12M,alpha1H,alpha2H,eiph,r,eipm1,eipm2,eipm12):  
    '''Main equations'''
    [SH, I1H, I2H, I12CH, I12SH, I1_2H, I2_1H, R1H, R2H, R12H,\
     SM, I1M, I2M, I12M, cumulative_I12CH, cumulative_I12SH, E1M, E2M,\
     E12M, E1H, E2H, E12CH, E12SH, E1_2H, E2_1H, I1E2M, I2E1M] = INP
    
    T = temp_func(tempMean,tempAmp,t)
    a = briere(T,13.35,40.08,2.02e-4)
    b1 = b2 = b12 = briere(T,17.05,35.83,8.49e-4)
    c1 = c2 = c12 = briere(T,12.22,37.46,4.91e-04)
    g = 1/(max(quadratic(T,9.16,37.73,-1.48e-01),0.01))
    m = gamma/g
 #   eipm = 1/(max(briere(T,10.68,45.9,6.65e-05),0.01))

    ##EIP values set here
 #   eipm1 = 2 #chikungunya
#    eipm2 = 14 #zika
#    eipm12 = 2 #chik, next will vary to be average of the two (8)

    lambda1H = m*a*b1*(I1M+I12M*p1H)
    lambda2H = m*a*b2*(I2M+I12M*p2H)
    lambda12H = m*a*b12*I12M*p12H
    
    lambda1M = a*c1*(I1H+I2_1H+(I12CH+I12SH)*p1M)
    lambda2M = a*c2*(I2H+I1_2H+(I12CH+I12SH)*p2M)
    lambda12M = a*c12*(I12CH+I12SH)*p12M
    
    dSH = -(lambda1H+lambda2H+lambda12H)*SH
    
    dE1H = lambda1H*SH - ((lambda2H + lambda12H)*(1-alpha1H) + (1./eiph))*E1H  
    dE2H = lambda2H*SH - ((lambda1H + lambda12H)*(1-alpha2H) + (1./eiph))*E2H  
    dE12CH = lambda12H*SH + lambda12H*((1-alpha2H)*E2H + (1-alpha1H)*E1H) - (1./eiph)*E12CH
    dE12SH = lambda2H*(1-alpha1H)*E1H + lambda1H*(1-alpha2H)*E2H - (1./eiph)*E12SH
    
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
    dSM = g*(E1M+E2M+E12M+I1M+I2M+I12M+I1E2M+I2E1M) - (lambda1M+lambda2M+lambda12M)*SM
    
    dE1M = lambda1M*SM - (lambda2M+lambda12M + g)*E1M - (1./eipm1)*E1M
    dE2M = lambda2M*SM - (lambda1M+lambda12M + g)*E2M - (1./eipm2)*E2M
    dE12M = lambda12M*SM + (lambda2M+lambda12M)*E1M + (lambda1M+lambda12M)*E2M - g*E12M - (1./eipm12)*E12M 

    dI1M = (1./eipm1)*E1M - lambda2M*I1M - g*I1M 
    dI2M = (1./eipm2)*E2M - lambda1M*I2M - g*I2M
    dI12M = (1./eipm12)*(E12M + I1E2M + I2E1M) - g*I12M
    dI1E2M = lambda2M*I1M - (1./eipm12)*I1E2M - g*I1E2M
    dI2E1M = lambda1M*I2M - (1./eipm12)*I2E1M - g*I2E1M
    
    dcumulative_I12CH = lambda12H*SH + lambda12H*((1-alpha2H)*I2H + (1-alpha1H)*I1H)
    dcumulative_I12SH = lambda2H*(1-alpha1H)*I1H + lambda1H*(1-alpha2H)*I2H
    
    Y = np.array([dSH, dI1H, dI2H, dI12CH, dI12SH, dI1_2H, dI2_1H, dR1H, dR2H,\
               dR12H, dSM, dI1M, dI2M, dI12M, dcumulative_I12CH, dcumulative_I12SH, dE1M, dE2M, \
               dE12M, dE1H, dE2H, dE12CH, dE12SH, dE1_2H, dE2_1H, dI1E2M, dI2E1M])
    
    Y[(Y<=0)*(INP<=0)]=0.0
    return Y   # For odeint
