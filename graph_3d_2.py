from __future__ import division
from src.GeneratedCarriers import LightGeneratedCarriers
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import pow

dt = 0.01
stepCnt = 30000

xs = np.empty((stepCnt + 1,))
a_s = np.empty((stepCnt + 1,))
b_s = np.empty((stepCnt + 1,))
c_s = np.empty((stepCnt + 1,))
d_s = np.empty((stepCnt + 1,))
e_s = np.empty((stepCnt + 1,))
f_s = np.empty((stepCnt + 1,))

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
xs[0], a_s[0], b_s[0], c_s[0], d_s[0], e_s[0], f_s[0]  = (V, 0., 0., 0., 0., 0., 0.)

for i in range(stepCnt):
    xs[i + 1] = xs[i] + 0.0001
    a_s[i + 1] = LightGeneratedCarriers(T, xs[i+1], 10**20/10**15, 10**20/10**15, 500*pow(10,-9), width, x, W, H_, intensity)
    b_s[i + 1] = LightGeneratedCarriers(T, xs[i+1], 10**20/10**15, 10**20/10**15, 550*pow(10,-9), width, x, W, H_, intensity)
    c_s[i + 1] = LightGeneratedCarriers(T, xs[i+1], 10**20/10**15, 10**20/10**15, 600*pow(10,-9), width, x, W, H_, intensity)
    d_s[i + 1] = LightGeneratedCarriers(T, xs[i+1], 10**20/10**15, 10**20/10**15, 650*pow(10,-9), width, x, W, H_, intensity)
    e_s[i + 1] = LightGeneratedCarriers(T, xs[i+1], 10**20/10**15, 10**20/10**15, 700*pow(10,-9), width, x, W, H_, intensity)
    f_s[i + 1] = LightGeneratedCarriers(T, xs[i+1], 10**20/10**15, 10**20/10**15, 750*pow(10,-9), width, x, W, H_, intensity)

fig = plt.figure()
#ax = fig.gca(projection='2d')
plt.plot(xs,a_s,'r--',xs,b_s,'b--',xs,c_s,'g--',xs,d_s,'y--',xs,e_s,'c--',xs,f_s,'k--')
# ax.set_xlabel("X Axis")
# ax.set_ylabel("Y Axis")
# ax.set_title("Doping")

plt.show()
