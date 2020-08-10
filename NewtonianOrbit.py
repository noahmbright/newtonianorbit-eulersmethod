import numpy as np
from  matplotlib import pyplot as plt


#Define parameters
G= 1.736e-29 #m^3/kg/yr^2
m1 = 5.972e24 #kg
m2 = 1.989e30 #kg
mu = m1*m2/(m1+m2) #kg
r0 = 1
phidot0 = np.pi*2  #1/yr
l = phidot0*mu*(r0**2)
gamma = G*m1*m2
params = gamma*mu/(l**2)



##define variables
t0 = 0
tf = 1
n=10000
deltat = (tf-t0)/(n-1)

t = np.linspace(t0, tf, n)
r = np.zeros([n])
v = np.zeros([n])

phi = np.zeros([n])


#Initial condition for dr/dt
v[0] = 0
#initial condition for r
r[0] = r0
#inital conditon for phi
phi[0] = 0


#Euler's method
for i in range(1,n):
    v[i] = v[i-1] + deltat*(l**2/mu/(r[i-1]**3) - gamma/(r[i-1]**2))/mu
    r [i] = r[i-1] + deltat*(v[i-1])
    phi[i] = phi[i-1] + deltat* (l/mu/r[i-1]**2)




graph = plt.subplot(111, projection = "polar")
graph.plot(phi,r)
graph.set_rmax=(2)
graph.set_rticks([0.5, 1])
graph.grid(True)
graph.set_title("Earth-Sun Separation Over 3 Years")
plt.show()
print(phi[len(phi)-1])