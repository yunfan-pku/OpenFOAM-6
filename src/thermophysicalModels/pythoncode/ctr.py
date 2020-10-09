#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 20:01:41 2020

@author: zhy
"""

import VLE
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
mix = VLE.solver(
    "/home/zhy/Documents/git/OpenFOAM-6/src/thermophysicalModels/testmain/system/thermo")
mix.specie = ["CH4", "O2"]
ch=0.001
mix.X = [ch, 1-ch]
mix.update()
delta = 0.2
dt=0.1
p = np.arange(40e5, 57e5, 0.1e5)
t = np.arange(140.0, 180.0, dt)

T, P = np.meshgrid(t, p)
V= T+3
Z = T+1
D = T+2
for ip in range(len(p)):
    mix.P=p[ip]
    for it in range(len(t)):
        mix.T = t[it]
        mix.solve(True)
        V[ip][it]=mix.vaporfra
        Z[ip][it]=mix.density()
for ip in range(len(p)):
    for it in range(1,len(t)):
        D[ip][it]=(Z[ip][it-1]-Z[ip][it])/dt


font1 = {'family' : 'Times New Roman',
'weight' : 'normal',
'size'   : 14,
}
cset0 =plt.contourf(t,p,V)
plt.colorbar(cset0)
plt.xlabel("$Temperature (K)$",font1)
plt.ylabel("Critical pressure (bar)",font1)
plt.show()
"""
cset =plt.contourf(t,p,Z)
plt.colorbar(cset)
plt.show()
cset1 =plt.contourf(t,p,D)
plt.colorbar(cset1)
plt.show()
"""