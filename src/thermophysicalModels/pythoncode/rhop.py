#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 20:57:06 2020

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
p = np.arange(45e5, 54e5, 1e5)
t = np.arange(140.0, 180.0, dt)
rho = np.arange(100, 800, 5.0)

Rho, P = np.meshgrid(rho, p)
V=P+2

TP= P+1
for ip in range(len(p)):
    mix.P=p[ip]
    for ir in range(len(rho)):
        TP[ip][ir]=mix.twophase(rho[ir],80.0,300.0)
        mix.solve(True)
        V[ip][ir]=mix.vaporfra


font1 = {'family' : 'Times New Roman',
'weight' : 'normal',
'size'   : 14,
}
cset =plt.contourf(rho,p,TP)
plt.xlabel("$Density (kg/m^3)$",font1)
plt.ylabel("pressure (bar)",font1)
plt.colorbar(cset)
plt.show()
"""
cset1 =plt.contourf(rho,p,V)
plt.colorbar(cset1)
plt.show()
"""
