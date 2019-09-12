from numpy import *
from scipy.integrate import odeint
from matplotlib.pyplot import *
import pandas as pd

Grdark=1./5000.

phi_data=array([0.22,1.8,6.2,13.5,18.0])*1e9
Gr_data=1./array([576.1,114.1,48.4,49.8,47.3])
phi0=1.8e8
ftphi=polyfit(log10(phi_data/phi0),Gr_data,1)
print ftphi
phi=logspace(8,10,100)


#Gr=Grdark+c*log(phi/phi0)

plot(log10(phi_data/phi0),Gr_data,'x')
plot(log10(phi_data/phi0),ftphi[0]*log10(phi_data/phi0)+ftphi[1],'r')
#semilogx(phi,Gr)
show()
