#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 20:57:06 2020

@author: zhy
"""


import VLE
import csv 
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
mix = VLE.solver(
    "/home/zhy/Documents/git/OpenFOAM-6/src/thermophysicalModels/testmain/system/thermo")
mix.specie = ["CH4", "O2"]
ch=0.3
mix.X = [ch, 1-ch]
mix.update()


ox=[]
op=[]
"""
for ix in np.arange(0., 1.1, 0.1):
    mix.X = [ix, 1-ix]
    irs=10.0
    for ip in np.arange(45e5,60e5,0.1e5):
        mix.P=ip
        pflag=False    
        for ir in np.arange(irs,800.0,1):
            twophaseflag=mix.twophase(ir,80.0,200.0)
            if twophaseflag==True:
                pflag=True
                irs=ir
                break
        if pflag==False:
            op.append(ip)
            ox.append(ix)
            print(ip,ix)
            break
"""
for ix in np.arange(0., 1.1, 0.1):
    mix.X = [ix, 1-ix]
    irs=10.0
    ip=40e5
    dp=10e5
    while dp>1:
        ip+=dp
        mix.P=ip
        pflag=False
        ir=irs
        while ir<800.0:
            twophaseflag=mix.twophase(ir,80.0,200.0)
            if twophaseflag==True:
                pflag=True
                irs=ir
                break
            ir+=5
        if pflag==False:
            ip-=dp
            dp/=2.0
    op.append(ip)
    ox.append(ix)
    print(ip,ix)

op2=np.array(op)
plt.plot(ox,op)

headers=["CH4","Pc"]
a=[ox,op]
b=[[row[i] for row in a] for i in range(len(a[0]))]
with open('/home/zhy/Documents/git/OpenFOAM-6/src/thermophysicalModels/pythoncode/test1.csv','w')as f:
    f_csv = csv.writer(f)
    f_csv.writerow(headers)
    f_csv.writerows(b)
        

