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
        #C1_d=prm['N']-O1-O2-C2
        #dC1dt=
        #prm=self.params
        #print(prm)
        """temporary zone
           replace with something more meaningful
           """
        F=prm['F'] #photons/ms
        pulse_start=10 #ms
        pulse_end=pulse_start+30#ms
        if t>pulse_start and t<pulse_end:
            Gr=1./5.#prm['Grdark']
        else:
            Gr=1./10000.#prm['Grdark']
        """end of temporary zone"""
        Ga1=prm['eps1']*F*self.f_singlepulse(t,prm['tau'],pulse_start,pulse_end)
        Ga2=prm['eps2']*F*self.f_singlepulse(t,prm['tau'],pulse_start,pulse_end)
        
        dO1dt=Ga1*C1-(prm['Gd1']+prm['e12'])*O1+prm['e21']*O2
        dO2dt=Ga2*C2-(prm['Gd2']+prm['e21'])*O2+prm['e12']*O1
        dC1dt=0#prm['Gd1']*O1+Gr*C2
        dC2dt=prm['Gd2']*O2-Ga2*C2-Gr*C2
        #dC1dt=-(dO1dt+dO2dt+dC2dt)
        
        dYdt=[dO1dt,dO2dt,dC1dt,dC2dt]
        return dYdt
    
    def __init__(self,csv_path='ChR2Params.csv'):
        params_tbl=pd.read_csv(csv_path)[['name','value']]
        par_keys=list(params_tbl['name'])
        par_vals=list(params_tbl['value'])
        self.params=dict(zip(par_keys,par_vals))
        

    def f_singlepulse(self,t,tau,tStart,tEnd):
        if t<=tStart:
            return 0
        elif t>tStart and t<tEnd:
            return 1-exp(-(t-tStart)/tau)
        elif t>=tEnd:
            return -exp( -(t-tStart)/tau)+exp(-(t-tEnd)/tau)



    
    def L2I(self,t,L):
        N=self.params['N']
        O1_0=0.0*N
        O2_0=0.0*N
        C1_0=0.61*N
        C2_0=N-O1_0-O2_0-C1_0
        Y0=[O1_0,O2_0,C1_0,C2_0]
        #print ("S1\t",Y0)
        #print ("S3\t",dYdt(t[0],Y0,self.params))
        Y=zeros([len(t),4])
        Y[0,:]=Y0
        for it,ct in enumerate(t[1:]):
            
            cY=list(odeint(self.dYdt,Y0,[ct,ct+self.params['dt']],args=(self.params,))[1].squeeze())
            #print ct,Y0,cY
            cY[2]=self.params['N']-cY[0]-cY[1]-cY[3]
            
            Y[it+1,:]=cY
            Y0=cY
        o1=Y[:,0]/self.params['N']
        o2=Y[:,1]/self.params['N']
        Imax=1 #nA
        I=self.params['Imax']*(o1+self.params['gamma']*o2)
        return I#odeint(self.dYdt,Y0,t,args=(self.params,))
        
ChR2=ChR2Dynamics()
tend=300

t=linspace(0,tend,tend*2+1)
for eps1 in [0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1]:
    ChR2.params['eps1']=eps1
    
    I=ChR2.L2I(t,100)
    plot(t+0.5,I)#,'k',linewidth=2)
#print (shape(Y))




#plot(t,o1,'.-')
#plot(t,o2,'.-')
#plot(t,c1,'.-')
#plot(t,c2,'.-')
#I=sum(Y,axis=1)

#legend(['100','10','1','0.1'])
#print (sum(Y,axis=1))
#print (ChR2.params)
show()
"""
t=arange(0,500,.2)
pulse_start,pulse_end=50,100
tau=30
f=array([f_singlepulse(ti,tau,pulse_start,pulse_end) for ti in t])
plot(t,f)
show()
"""
