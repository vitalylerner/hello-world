# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 11:28:43 2019

@author: Vitaly Lerner
"""
from numpy import *
from scipy.integrate import odeint
from matplotlib.pyplot import *
import pandas as pd




class ChR2Dynamics:
    
    params={}
    
    def dYdt(self,Y,t,prm):
        """Source:Nikolic et al, Photocycles of channelrhodopsin-2"""
        #print ("dYdt t\t" , t)
        #print ("dYdt Y\t" , Y)
        #print ("dYdt 1\t" , prm)
        O1,O2,C1,C2=Y
        """if O1>prm['N']:
            O1=prm['N']
        if O2>prm['N']:
            O2=prm['N']
        if C1>prm['N']:
            C1=prm['N']
        if C2>prm['N']:
            C2=prm['N']"""
        C1=prm['N']-O1-O2-C2
        #prm=self.params
        #print(prm)
        """temporary zone
           replace with something more meaningful
           """
        F=prm['F'] #photons/ms
        pulse_start=100 #ms
        pulse_end=pulse_start+200#ms
        if t>pulse_start and t<pulse_end:
            Gr=1./5.#prm['Grdark']
        else:
            Gr=1./10000.#prm['Grdark']
        """end of temporary zone"""
        Ga1=prm['eps1']*F*f_singlepulse(t,prm['tau'],pulse_start,pulse_end)
        Ga2=prm['eps2']*F*f_singlepulse(t,prm['tau'],pulse_start,pulse_end)
        
        dO1dt=Ga1*C1-(prm['Gd1']+prm['e12'])*O1+prm['e21']*O2
        dO2dt=Ga2*C2-(prm['Gd2']+prm['e21'])*O2+prm['e12']*O1
        dC1dt=0#prm['Gd1']*O1+Gr*C2
        dC2dt=prm['Gd2']*O2-Ga2*C2-Gr*C2
        
        dYdt=[dO1dt,dO2dt,dC1dt,dC2dt]
        return dYdt
    
    def __init__(self,csv_path='ChR2Params.csv'):
        params_tbl=pd.read_csv(csv_path)[['name','value']]
        par_keys=list(params_tbl['name'])
        par_vals=list(params_tbl['value'])
        self.params=dict(zip(par_keys,par_vals))
        

    def f_singlepulse(t,tau,tStart,tEnd):
        if t<tStart:
            return 0
        elif t>=tStart and t<tEnd:
            return 1-exp(-(t-tStart)/tau)
        elif t>=tEnd:
            return -exp( -(t-tStart)/tau)+exp(-(t-tEnd)/tau)



    
    def L2I(self,t,L):
        N=self.params['N']
        O1_0=0
        O2_0=0
        C1_0=N/1000
        C2_0=N-O1_0-O2_0-C1_0
        Y0=[O1_0,O2_0,C1_0,C2_0]
        print ("S1\t",Y0)
        #print ("S3\t",dYdt(t[0],Y0,self.params))
        return odeint(self.dYdt,Y0,t,args=(self.params,))
        
ChR2=ChR2Dynamics()
t=linspace(0,1000,1001)

for F in [100,10,1,0.1]:#linspace(1,,50):
    ChR2.params['F']=F
    Y=ChR2.L2I(t,100)
    
    #print (shape(Y))
    gamma=0.05
    o1=Y[:,0]/ChR2.params['N']
    o2=Y[:,1]/ChR2.params['N']#,C1,C2=Y
    Imax=1 #nA
    I=Imax*(o1+gamma*o2)
    
    #plot(t,o1,'.-')
    #plot(t,o2,'.-')
    
    plot(t,I)
legend(['100','10','1','0.1'])
#print (sum(Y,axis=1))
#print (ChR2.params)

"""
t=arange(0,500,.2)
pulse_start,pulse_end=50,100
tau=30
f=array([f_singlepulse(ti,tau,pulse_start,pulse_end) for ti in t])
plot(t,f)
show()
"""