#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 17:27:13 2020

@author: zhy
"""

import csv 
import matplotlib.pyplot as plt
data=list(csv.reader(open("test.csv")))[1:]
font1 = {'family' : 'Times New Roman',
'weight' : 'normal',
'size'   : 14,
}
x=[]
p=[]
for i in range(len(data)):
    x.append(float(data[i][0]))
    p.append(float(data[i][1])/1e5)
plt.plot(x,p)
plt.xlabel("$x_{CH_4}$",font1)
plt.ylabel("Critical pressure (bar)",font1)
plt.savefig("CH4_O2.png",dpi=300)