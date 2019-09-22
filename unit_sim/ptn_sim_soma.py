# -*- coding: utf-8 -*-
"""
Vitaly Lerner, 2019
Simulation of a single neuron expressing ChR2
to a pattern of light activations

Functions for numerical integration of light on soma
"""
from numpy import *

def sphere_equi(n=500):
    #distribute points represented by Eucledean coordinates
    #use deterministic method (non-random)
    #achieve equidistribution with the method described there:
    ##https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
    a=4*pi/n
    d=sqrt(a)
    Mth=int(round(pi/d))
    #print Mth
    dth=pi/Mth
    dphi=a/dth

    k=0
    theta=[]
    phi=[]
    for m in range(Mth):
        cTh=pi*(m+0.5)/Mth
        Mphi=int(round(2*pi*sin(cTh)/dphi))
        for n in range(Mphi):
            cPhi=2*pi*n/Mphi
            theta.append(cTh)
            phi.append(cPhi)
            

    x = sin(theta) * cos(phi)
    y = sin(theta) * sin(phi)
    z = cos(theta)
    return x,y,z
 
    
def soma_integrate_light(r_soma,r_spot,x_spot):
    v = exp( - ((r_soma*crd_x-x_spot)**2+(r_soma*crd_y)**2)/(2*r_spot**2) )
    return sum(v)*4*pi*r_soma**2/len(v)