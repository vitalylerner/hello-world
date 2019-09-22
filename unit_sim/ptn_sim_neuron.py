# -*- coding: utf-8 -*-
"""
Vitaly Lerner, 2019
Simulation of a single neuron expressing ChR2
to a pattern of light activations
Geometry: the slices are coronal, from the right hemisphere
    the coordinates are experiment-oriented

    z: rostral-caudal, z=0 is the top surface of the slice,
        i.e. nose --> 0um ---300um -->  tail
    y: dorsal-ventral, y=0 is the pia, i.e.
        dura --> pia=0um --->900um --> white matter
    x: arbitrary coordinate orthogonal to (y,z) such that approsimately
       x axis is parallel to a local line of pia and wm
       approximately right to left
"""
import pandas as pd
from numpy import *
from matplotlib.pyplot import *
import os,shutil,pickle
from multiprocessing import Pool,Process
import time
from mpl_toolkits.mplot3d import axes3d, Axes3D
from scipy import integrate
from ChR2Dynamics import ChR2Dynamics
from ptn_sim_morphology import *
from ptn_sim_optostim import optostim
from ptn_sim_soma import sphere_equi,soma_integrate_light


SEGMENT_LENGTH=8 #micron
swc_path_base_allen='../AllenCells/SWC/'
def list_local_cells_allen():
    return [int(fn[:-4]) for fn in os.listdir(swc_path_base_allen)]
    


n_points_on_sphere=500 #keep constant
crd_x,crd_y,crd_z=sphere_equi(n_points_on_sphere)


"""
End of Light on Soma integration
"""

class neuron:
    # parameters are stored as pandas dataframe
    # of structure
    # nams, datatype, units, value
    
    """procedure for importing swc and exporting a cell
    n.assign_Allen_ID()/n.assign_Lerner_ID()
    n.assign_swc(path_to_swc)
    n.export_morphology()
    """
    params=None        #morphological and physiological parameters
    neuron_id=None    #identifier
    morph=None        #morphology
    
    
    def export_morphology(self):
        src=self.get_param_value('morph_source')
        cellid=self.get_param_value('morph_cell_id')
        
        base_name=src+'_{}'.format(cellid)
        fPickle='MORPH/'+base_name+'/morph.pkl'
        
        self.morph.export_csv(base_name)
        pickle.dump(self.morph,open(fPickle, "wb" ))

    def import_morphology(self):
        src=self.get_param_value('morph_source')
        cellid=self.get_param_value('morph_cell_id')
        
        base_name=src+'_{}'.format(cellid)
        self.morph.import_csv(base_name)
        
        
    def get_param_units(self,name):
    #retrieve units of a certain parameter
        sp=self.params
        cU=sp[sp['name']==name]['units'].to_string(index=False)
        return cU

    def get_param_value(self,name):
    #retrieve a value of a certain parameter
    #input: parameter name, string
    #output: parameter value, parameter-dependent datatype
        sp=self.params
        cValue=sp[sp['name']==name]['value']
        cDType=sp[sp['name']==name]['datatype'].to_string(index=False)
        cValue2=cValue.to_string(index=False)
        if cDType=='bool':
            cVal=cValue2=='True'
        elif cDType=='real':
            cVal=float(cValue2)
        elif cDType=='string':
            cVal=cValue2
        elif cDType=='int':
            cVal=int(cValue2)
        return cVal
        
    def set_param_value(self,name,value):
        i=self.params[self.params['name']==name].index
        self.params.ix[i,'value']=value

    def assign_morphology(self,M):
        self.morph=M
        
    def assign_swc(self,swc_path):
        M=morphology(swc_path)
        self.assign_morphology(M)
        
    def assign_Allen_ID(self,CellID_Allen):
        self.set_param_value('morph_source','ALLEN')
        self.set_param_value('morph_cell_id',CellID_Allen)
        #print self.params
        
    def assign_Lerner_ID(self,CellID_Lerner):
        self.set_param_value('morph_source','LERNER')
        self.set_param_value('morph_cell_id',CellID_Lerner)
        #print self.params
        
    def list_params(self):
    #retrieve a list of all physiological and morphological parameters
    #output: list of strings
        sp=self.params
        cColNames=list(sp[:]['name'])
        return cColNames

    def __init__(self,neuron_id,params_csv='neuron_default.csv'):
    #initiator
    #neuron_id: identificator
    #params_csv: path to a list of all physiological and morphologilcal
    #            parameters
        dt = pd.read_csv(params_csv)
        self.params=dt
        self.neuron_id=neuron_id
        self.morph=morphology()
    
    def spot_flux(self,spot):
        """
        Apply a single spot (timeless) to a neuron
        Input: spot : dict of x,y,r
        Output: table of segments: dist (um), light flux (au)
        """
        tblSegs=self.morph.segments_geometry()
        l=array(tblSegs['length'])
        d=array(tblSegs['dist'])
        
        sCoM,rSoma=self.morph.soma_geometry()
        xSoma=sCoM[0]
        ySoma=sCoM[1]
        
        x0=array(tblSegs['x0'])
        y0=array(tblSegs['y0'])
        
        x1=array(tblSegs['x1'])-x0
        y1=array(tblSegs['y1'])-y0
        
        xs=spot['x']-x0
        ys=spot['y']-y0
        rs=spot['r']
        x0t=0*x0
        y0t=0*y0
        
        dx=x1-x0
        dy=y1-y0
        k=dy/dx
        alpha=arctan(k)
        
        x0t=x0*cos(alpha)+y0*sin(alpha)
        y0t=-x0*sin(alpha)+y0*cos(alpha)
        
        x1t=x1*cos(alpha)+y1*sin(alpha)
        y1t=-x1*sin(alpha)+y1*cos(alpha)
        
        xst=xs*cos(alpha)+ys*sin(alpha)
        yst=-xs*sin(alpha)+ys*cos(alpha)
        
        rNeg=x1t<0
        x1t[rNeg]=-x1t[rNeg]
        xst[rNeg]=-xst[rNeg]
        yst[rNeg]=-yst[rNeg]
        

        L=lambda x,y,xs,ys,rs:1.*exp(- ((x-xs)**2+(y-ys)**2)/(2*rs**2))
        #ax=gca()
        #ax.add_artist(Circle((spot['x'],spot['y']),spot['r'],alpha=0.6,color=[0.2,0.2,0.8]))
        NSeg=len(x0)
        Flux=zeros(NSeg)
        for iSeg in range(NSeg):#[10,13,14]:
            flux=integrate.dblquad(L, 0, x1t[iSeg], lambda x: -.1*l[iSeg], lambda x: .1*l[iSeg],args=(xst[iSeg],yst[iSeg],rs))[0]+0.001
            Flux[iSeg]=flux
            #clrval=0.05+1.2*(log10(flux)+3)/5
            #plot( [x0[iSeg],x1[iSeg]+x0[iSeg]],[y0[iSeg],y1[iSeg]+y0[iSeg]],color=clrval*ones(3))
            
        #separately calculate soma
        dx_soma=spot['x']-xSoma
        dy_soma=spot['y']-ySoma
        
        delta_soma=sqrt(dx_soma**2+dy_soma**2)
        
        return pd.DataFrame(data={'flux':Flux,'dist':d})
        
    def apply_optostim(self,lstim):
        #Parameters of light stimuli
        dur=lstim.seq_params['dur']
        ISI=lstim.seq_params['ISI']
        lst=lstim.seq_params['lst']
        nstim=len(lst)
        ChR2=ChR2Dynamics()
        dur_tot=nstim*ISI+dur+200
        
        r=lstim.spots_r
        
        #morphology
        Nseg=len(array(self.morph.segments['segment_id']))
        
        #Light generation
        L=zeros([Nseg,dur_tot])
        
        
        #superimpose light flux from all stimuli spots
        for isp,nsp in enumerate(lst[:250]):
             cspot={'x':lstim.spots_xy[nsp,0],'y':lstim.spots_xy[nsp,1],'r':r}
             ct=isp*ISI
             fld=self.spot_flux(cspot)
             flux=fld['flux']
             dist=fld['dist']
             tmp=outer(ones([1,dur]),flux).T
             #print 'A',shape(tmp),shape( L[:,ct:ct+dur])
             L[:,ct:ct+dur]+=tmp
             
        #translate light flux to current for each segment separately
        dt=ChR2.params['dt']
        nrep=int(1./dt)
        L=repeat(L,nrep,axis=1)
        I=0*L
        t=arange(shape(L)[1])*dt
        for iSeg in range(Nseg):
            l=L[iSeg,:].squeeze()
            I[iSeg,:]=ChR2.L2I(t,l)
        return I+10
        #print L
        #matshow(log10(L))
        #show()



       
lstCells=list_local_cells_allen()
N=[None for id in lstCells]
print (time.asctime( time.localtime(time.time()) ))
for iN,CellID in enumerate(lstCells[:1]):
    N[iN]=neuron(iN)
    fSWC=swc_path_base_allen+'{}.swc'.format(CellID)
    N[iN].assign_Allen_ID(CellID)
    N[iN].import_morphology()
    N[iN].morph.translate([15*iN,0,0])
    #N[iN].morph.draw({'Layout':'Branches','alpha':0.4})
    #N[iN].morph.draw({'Layout':'Soma','alpha':0.4})
    #N[iN].spot_flux({'x':-50,'y':-600,'r':20})
    
    #N[iN].morph.soma_geometry()

print (time.asctime( time.localtime(time.time()) ))

map_params={'r':10,'dr':20,'x0':0,'y0':-570,'N':3}
lstim=optostim(map_params)
S=[8,18,16,13,10]
lstim.assign_seq({'lst':S,'ISI':5,'dur':12})

I=N[0].apply_optostim(lstim)

print (time.asctime( time.localtime(time.time()) ))
#N[0].morp.draw({
#N[0].apply_spot({'x':-50,'y':-700,'r':20})
IM=log10(I+0.1)
#IM[IM<0]=0
#IM=IM/IM.max()
matshow(IM,cmap='hot')
print (time.asctime( time.localtime(time.time()) ))
#axis('equal')
show()
