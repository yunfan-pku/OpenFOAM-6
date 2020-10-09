#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 21:09:24 2020

@author: zhy
"""
import VLE
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
 
delta = 0.2
p = np.arange(40e5, 57e5, 0.01e5)
t = np.arange(100.0, 180.0, 0.1)
tl=100
tr=180
pl=1e5
pr=180e5
T, P = np.meshgrid(t, p)
Z = T**2 + P**2

mix = VLE.solver(
    "/home/zhy/Documents/git/OpenFOAM-6/src/thermophysicalModels/testmain/system/thermo")
mix.specie = ["CH4", "O2"]
ch=0.001
mix.X = [ch, 1-ch]
mix.P = 500e5
mix.update()
mix.solve(True)
cp=[]
cx=[]
tcp=0.0
for ix in np.arange(0.1, 0.994, 0.03):
    mix.X = [ix, 1-ix]
    tcp=0.0
    pl_temp=pl
    pr_temp=pr
    while pr_temp-pl_temp>0.01e5:
            pflag=False
            pm=(pl_temp+pr_temp)/2
            mix.P=pm
            tl_temp=tl
            tr_temp=tr
            while tr_temp-tl_temp>3e-4:
                tm=(tl_temp+tr_temp)/2
                mix.T = tm
                mix.solve(True)
                if(mix.vaporfra*(1-mix.vaporfra)>0.04):
                    if(pm>tcp):
                        tcp=pm
                        pflag=True
                        break
                else:
                    if(mix.vaporfra>0.5):
                        tr_temp=tm
                    else:
                        tl_temp=tm
            if(pflag==True):
                pl_temp=pm
            else:
                pr_temp=pm
                    
    if(tcp>0):
        cx.append(ix)
        cp.append(tcp)
            
            #Z[ip][it]=mix.vaporfra*(1-mix.vaporfra)#mix.density()

plt.plot(cx,cp)          
"""
cset =plt.contourf(t,p,Z)
plt.colorbar(cset)
"""
"""
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(t, p, z, cmap=cm.jet, linewidth=0.01)
plt.show()
"""