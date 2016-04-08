from __future__ import division
from src.GeneratedCarriers import LightGeneratedCarriers, LengthElectronDiffusion, LengthHoleDiffusion, Area
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import pow, exp

dt = 0.0001
stepCnt = 20000

def formula(T, V, mn, mp, xp, xn, d):
    A = 1
    B = 1
    C = 1
    D = 1
    f = 0
    integr1 = (10**20/10**15*(-d/2+xp) + D*exp((d/2+xp)/LengthElectronDiffusion(T))*LengthElectronDiffusion(T) - D*LengthElectronDiffusion(T))
    temp1 =  C*(LengthElectronDiffusion(T)*exp((-1)*(d/2-xp)/LengthElectronDiffusion(T))) - C*LengthElectronDiffusion(T)
    integr1  = integr1 + temp1

    integr2 = (10**20/10**15*(d/2-xp) + B*exp((d/2-xn)/LengthHoleDiffusion(T))*LengthHoleDiffusion(T)) - B*LengthHoleDiffusion(T)
    temp2 =  A*(LengthHoleDiffusion(T)*exp((-1)*(x-xn)/LengthHoleDiffusion(T))) - A*LengthElectronDiffusion(T)
    integr2  = integr2 + temp2

    return ((1.602*(10**(-19))*Area(width))/(d**2))*(f-V)*(mp*integr2+mn*integr1)

xs = np.empty((stepCnt + 1,))
a_s = np.empty((stepCnt + 1,))

# Setting initial values
V = -1.0000
ni = 10**10
nd = 10**15
T = 300
wavelength = 800*pow(10,-9)
width = pow(10,-3)
W = pow(10,-6)
x = 2*pow(10,-3)
H_ = 5*pow(10,-3) - W/2
intensity = 1345
mn = 0.14
mp = 0.045
xn = 0.613*(10**(-6))
xp = 0.613*(10**(-6))
d = 0.001;
xs[0], a_s[0] = (V, formula(T,V,mn,mp,xp,xn,d))

for i in range(stepCnt):
    xs[i + 1] = xs[i] + 0.0001
    a_s[i + 1] = (10**(-18))*(exp(xs[i+1]/0.7)-1)-formula(T,xs[i+1],mn,mp,xp,xn,d)

fig = plt.figure()
#ax = fig.gca(projection='2d')
plt.plot(xs,a_s,'r--')
# ax.set_xlabel("X Axis")
# ax.set_ylabel("Y Axis")
# ax.set_title("Doping")

plt.show()
