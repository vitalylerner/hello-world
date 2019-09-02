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
        
n1=neuron(0)
#a=n1.get_param('tau')
c=n1.list_params()
for p in c: 
	a=n1.get_param_value(p)
	u=n1.get_param_units(p)
	print (p,'\t',a,'\t',u)

