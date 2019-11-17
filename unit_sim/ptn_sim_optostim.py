# -*- coding: utf-8 -*-
"""
Vitaly Lerner, 2019
Simulation of a single neuron expressing ChR2
to a pattern of light activations
Optostim: Optical stimulation in 2D
          produced using DMD 
    

"""

import pandas as pd
from numpy import *
from matplotlib.pyplot import *

class optostim:
    spots_xy=[]
    spots_r=0
    map_params={}
    seq_params={}
    
    def circ_nspots(self,n):
        if n>5:
            return 0
        else:
            return [1,8,12,18,25,31][n]
            
    def __init__    (self,map_params=None):
        if not map_params==None:
            self.build_map(map_params)
    
    def assign_seq(self,seq_params):
        self.seq_params=seq_params
        
    
    def draw        (self,list_special=[]):
        Map=self.spots_xy
        r=self.spots_r
        for im in range(shape(Map)[0]):
            m=Map[im,:]
            plot(m[0],m[1],'.',alpha=0)
            ax=gca()
            ax.add_artist(Circle((m[0],m[1]),r,alpha=0.1,color=[0.2,0.2,0.8],linewidth=0))
            if im in list_special:
                text(m[0],m[1],'{}'.format(im),horizontalalignment='center',   verticalalignment='center',size=9,alpha=0.7,color=[0.9,0.2,0.2])
            else:
                text(m[0],m[1],'{}'.format(im),horizontalalignment='center',   verticalalignment='center',size=8,alpha=0.5)
        
    def build_map   (self,map_params):
        self.map_params=map_params
        r=map_params['r']
        dr=map_params['dr']
        x0=map_params['x0']
        y0=map_params['y0']
        N=map_params['N']
        theta0=map_params['theta']
        self.spots_r=r
        spot0=[x0,y0]
        Map=array(spot0)
        
        for iCircle in range(1,N+1):
            c_n=self.circ_nspots(iCircle)
            c_r=dr*iCircle
            c_dth=2*pi/c_n
            c_th=arange(c_n)*c_dth+theta0
            c_x=c_r*cos(c_th)+x0
            c_y=c_r*sin(c_th)+y0
            cp=vstack([c_x,c_y]).T
            Map=vstack([Map,cp])
        self.spots_xy=Map
