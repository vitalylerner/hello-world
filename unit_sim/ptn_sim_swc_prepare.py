import pandas as pd
from numpy import *
from matplotlib.pyplot import *
import os,shutil,pickle,sys

if sys.version_info[0]<3:
    PYTHON=2
else:
    PYTHON=3
if os.name=='nt':
    WINDOWS=True
    #print 'WINDOWS!!'
else:
    WINDOWS=False


swc_folder='../MyCells/20190724_2/'
swc_export='../MyCells/SWC/201907242.swc'

swc_files=sort(os.listdir(swc_folder))
iS=0
D=[]
for iSWC,fSWC in enumerate(swc_files):
    swc_dt=pd.read_csv(swc_folder+fSWC, skiprows=6,sep='\s',header=None)
    swc_dt.columns=['id','type','x','y','z','r','pid']
    swc_dt['id']=swc_dt['id']+iS
    swc_dt['pid']=swc_dt['pid']+iS*(swc_dt['pid']>=0)+2*(swc_dt['pid']==-1)
    if iSWC==0:
        swc_dt['type']=4
    else:
        swc_dt['type']=3
    D.append(swc_dt)
    #swc_dt[swc_dt['pid']==0]
    iS+=len(swc_dt)
    
swc_header = [
         '# Traced by Vitaly Lerner\n# using simple neurite tracer in FIJI and united later using python\n# id type x y z r pid',
        '','','','','','']
#print swc_header

swc=pd.concat(D,ignore_index=True)
swc.loc[0,'type']=1
swc[['x','y']]=swc[['x','y']]/13*520*1.154/5.*6.
swc.to_csv(swc_export,index=False,sep=' ',header=swc_header)

    
        
