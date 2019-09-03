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


class neuron:
	# parameters are stored as pandas dataframe
	# of structure 
	# nams, datatype, units, value
	
	params=None		#morphological and physiological parameters
	neuron_id=None	#identifier
	
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
		if cDType=='bool':
			cValue2=cValue.to_string(index=False)
			cVal=cValue2=='True'
		elif cDType=='real':
			cValue2=cValue.to_string(index=False)
			cVal=float(cValue2)
		return cVal
    
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

class morphology:
	swc=None
	swc0=None
	
	def read_swc (self,fSWC):
		#read standard morphology ascii file and 
		#return: pandas dataframe of swc
		swc_dt=pd.read_csv(fSWC, skiprows=3,sep='\s',header=None)
		swc_dt.columns=['id','type','x','y','z','r','pid']
		return swc_dt
		
	def assign(self, fSWC):
		#assign morphology from file
		#store a copy of orginal swc data for reference
		self.swc=self.read_swc(fSWC)
		self.swc0=self.read_swc(fSWC)
		
	def align(self):
		#move the neuron to (0,0)
		#rotate it such that its apical dendrite would reach uniform pia
		#works for pyramidal neurons,
		#won't work for PV
		dt=self.swc
		#extract x,y coordinates of apical dendtite
		xy0=array(dt[:][['x','y']])
		x0=xy0[:,0]
		y0=xy0[:,1]
		
		xy_apical=array(dt[dt['type']==4][['x','y']])
		xy_soma=array(dt[dt['type']==1][['x','y']])

		#extract x,y coordinates and place soma at (0,0)
		x = xy_apical[:,0]-xy_soma[0,0]
		y = xy_apical[:,1]-xy_soma[0,1]

		#Rough estimation of rotation angle
		p1=polyfit(x,y,1)
		theta=-arctan(p1[0])
		x1=x*cos(theta)-y*sin(theta)
		y1=x*sin(theta)+y*cos(theta)
		if mean(x1)<0:
			x1=-x1
			y1=-y1
			theta+=pi
			
		#fine estimation of rotation angle
		rng_middle=[i for i,cx in enumerate(x1) if cx>150 and cx<450]
		x2=x1[rng_middle]
		y2=y1[rng_middle]
		p2=polyfit(x2,y2,1)
		theta2=-arctan(p2[0])

		Theta=theta+theta2+pi/2
		x3=x0*cos(Theta)-y0*sin(Theta)
		y3=x0*sin(Theta)+y0*cos(Theta)
		y3-=mean(sorted(y3)[-20:])
		x3-=x3[0]
		self.swc[:]['x']=x3
		self.swc[:]['y']=y3
		
	def __init__(self,fSWC=None,bAlign=True):
		#Constructor
		if not fSWC==None:
			self.assign(fSWC)
			if bAlign:
				self.align()
			
	def translate (self, v):
		#move by a (x,y,z) vector
		self.swc[:]['x']+=v[0]
		self.swc[:]['y']+=v[1]
		self.swc[:]['z']+=v[2]
		
	def draw(self,params=None):
		#draw using pyplot library
		if params==None:
			style='.k'
		else:
			style=params['style']
		xy0=array(self.swc[:][['x','y']])
		x0=xy0[:,0]
		y0=xy0[:,1]
		plot(x0,y0,style,markersize=.2)
		axis('equal')
		
#print dt.head()

F=['471129934.swc','515524026.swc','515249852.swc']
S=[{'style':'.'},{'style':'.'},{'style':'.'}]
#fSWC='471129934.swc'
#fSWC='515524026.swc'
#fSWC='515249852.swc'
for fSWC,dp,iN in zip(F,S,range(3)):
	M=morphology(fSWC)
	
	M.align()
	M.translate((100*iN,0,0))
	M.draw(dp)
show()
"""
n1=neuron(5)
#a=n1.get_param('tau')
c=n1.list_params()
for p in c: 
	a=n1.get_param_value(p)
	u=n1.get_param_units(p)
	print (p,'\t',a,'\t',u)

"""
