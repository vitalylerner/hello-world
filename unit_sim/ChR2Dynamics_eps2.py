from numpy import *
from scipy.integrate import odeint
from matplotlib.pyplot import *
import pandas as pd

Grdark=1./5000.

phi_data=array([0.18,1.8,18.0])*1e9
eps2_data=1./array([0.01,0.03,0.1])
phi0=1.75e8


ftphi=polyfit(log10(phi_data/phi0),log10(eps2_data),1)
print ftphi
#phi=logspace(8,10,100)


#Gr=Grdark+c*log(phi/phi0)

plot(log10(phi_data/phi0),log10(eps2_data),'x-')
plot(log10(phi_data/phi0),ftphi[0]*log10(phi_data/phi0)+ftphi[1],'r')
#semilogx(phi,Gr)
show()
