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
    def f(self,t):
        f_t=self.params['ft_t']
        f_f=self.params['ft_f']
        
        return f_f[argmin(abs(f_t-t))]
        
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
    def L2f(self,t,L,tau_ChR2=1.3):
        #Convert light flux arbitrary function of 
        #the form of steps into 
        #the f(t,t_ChR2) 
        E=list(where(hstack([[False],abs(L[1:]-L[:-1])>0]))[0])
        F=t*0
        #plot(t,L)
        tau=1.3
        Eend=E[1:]
        Eend.append(len(t)-1)
        for indEdge,indEdgeNext in zip(E,Eend):
            l0=L[indEdge-1]
            l1=L[indEdge+1]
            t0=t[indEdge]
            t1=t[indEdgeNext]
            if l1>l0:
                cf=(l1-l0)*(1-exp(- (t-t0)/tau ))+l0
            else:
                cf=(l0-l1)*exp (-(t-t0)/tau)+l1
            F[indEdge:indEdgeNext]=cf[indEdge:indEdgeNext]
        return F
        
    def dYdt(self,Y,t,prm,fval,eps2,Gr):
        """Source:Nikolic et al, Photocycles of channelrhodopsin-2"""
        #print ("dYdt t\t" , t)
        #print ("dYdt Y\t" , Y)
        #print ("dYdt 1\t" , prm)
        O1,O2,C1,C2=Y
        """temporary zone
           replace with something more meaningful
           """
        F=prm['F'] #photons/ms
        Ga1=prm['eps1']*F*fval
        Ga2=eps2*F*fval
        
        dO1dt=Ga1*C1-(prm['Gd1']+prm['e12'])*O1+prm['e21']*O2
        dO2dt=Ga2*C2-(prm['Gd2']+prm['e21'])*O2+prm['e12']*O1
        dC1dt=0
        dC2dt=prm['Gd2']*O2-Ga2*C2-Gr*C2

        
        dYdt=[dO1dt,dO2dt,dC1dt,dC2dt]
        return dYdt


    def L2I(self,t,L):
        prm=self.params
        #initial conditions
        N=prm['N']
        O1_0=0.0*N
        O2_0=0.0*N
        C1_0=0.61*N
        C2_0=N-O1_0-O2_0-C1_0
        Y0=[O1_0,O2_0,C1_0,C2_0]

        Y=zeros([len(t),4])
        Y[0,:]=Y0
        #Calculate light-dynamics-dependent functions
        
        #f(t,tau_ChR2)
        ft=self.L2f(t,L,self.params['tau'])
        
        #Gr(phi)
        phi=L*prm['phimax']/prm['phi0']
        phi_zero=(phi==0)

        Gr=prm['Grdark']+prm['Grslope']*log10(phi)
        Gr[phi_zero]=prm['Grdark']
        
        #eps2(phi)
        eps2=prm['eps2_p0']*log10(phi)+prm['eps2_p1']
        eps2[eps2>0.5]=.5
        eps2[eps2<0.001]=.001

        for it,ct in enumerate(t[1:]):
            
            
            cY=list(odeint(self.dYdt,Y0,
                    [ct,ct+self.params['dt']],
                    args=(self.params,ft[it],eps2[it],Gr[it]))[1].squeeze())
            
            #Apply boundary condition of preserving
            #the number of channel copies
            cY[2]=self.params['N']-cY[0]-cY[1]-cY[3]
            
            Y[it+1,:]=cY
            
            Y0=cY
        o1=Y[:,0]/self.params['N']
        o2=Y[:,1]/self.params['N']
        I=self.params['Imax']*(o1+self.params['gamma']*o2)
        return I
       
"""TESTING SECTION
ChR2=ChR2Dynamics()
tend=200

def stepf(t,tstart,tend):
    if hasattr(t,'__len__'):
        return array([ct>tstart and ct<=tend for ct in t])
    elif isinstance(t,float):
        return float(t>tstart and t<=tend)
    else:
        return None

def dummy_light(tend=200,dt=0.5,dur=40):
    ps=[10,30,50,70]
    pA=[0.1,0.3,0.5,0.7]#[::-1]
    t=linspace(0,tend,tend/dt+1)
    L=t*0
    for ips,cps in enumerate(ps):
        L+=pA[ips]*stepf(t,cps,cps+dur)
    return t,L
"""

