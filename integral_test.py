from matplotlib.pyplot import *
from numpy import *
from scipy import integrate
from mpl_toolkits.mplot3d import axes3d, Axes3D
from mpl_toolkits import mplot3d
import time
#import random





def soma_integrate_light(r_soma,r_spot,x_spot):

    v = exp( - ((r_soma*crd_x-x_spot)**2+(r_soma*crd_y)**2)/(2*r_spot**2) )

    return sum(v)*4*pi*r_soma**2/len(v)



n=500
r_soma=10
r_spot=20
x_spot=50

X_spot=arange(-100,100,5)
R_spot=arange(0.1,200,0.01)

#print X_spot


#I=[soma_integrate_light(r_soma,r_spot,x_spot) for x_spot in X_spot]
print time.asctime( time.localtime(time.time()) )
I=[soma_integrate_light(r_soma,r_spot,x_spot) for r_spot in R_spot]
print time.asctime( time.localtime(time.time()) )
print len(I)
#plot(X_spot,I)
plot(R_spot,I)
show()
"""
#print x
#    
#f = lambda y, x,a1: a1
#x=integrate.dblquad(f, 0, 1, lambda x: 0, lambda x: 1,args=[5])[0]
#print x

f=figure()
#ax=Axes3D(f)

clr=v-min(v)
clr=clr/max(clr)
for cx,cy,cz,cc in zip(x,y,z,clr):
    #ax.scatter(cx,cy,cz,'.',color=cc*ones(3))
    plot(cx,cy,'.',color=cc*ones(3))
#plot(theta,phi,'.')
xlabel('x')
show()
"""


