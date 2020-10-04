# -*- coding: utf-8 -*-
"""
Spyder 编辑器

这是一个临时脚本文件。
"""
import VLE
import numpy as np
import matplotlib.pyplot as plt
import csv 
mix = VLE.solver(
    "./src/thermophysicalModels/testmain/system/thermo")
mix.specie = ["CH4", "O2"]
mix.X = [0.4, 0.6]
mix.T = 250
mix.P = 80e5
mix.update()
mix.solve(True)
p=[]
x=[]
tp=40e5;
lp=0;
for ix in np.arange(0.01, 0.99, 0.01):
    mix.X = [float(ix), float(1.0-ix)]
    lp=tp
    tp=0
    for ip in np.arange(lp-0.5e5, 53e5, 0.01e5):
        mix.P=float(ip)
        flag=0
        for i in np.arange(80, 200, 0.1):
            mix.T = float(i)
            mix.solve(True)
            if(mix.vaporfra < 0.999999 and mix.vaporfra>2e-7):
                temp=mix.equalconstant
                #print(i)
                #print(mix.vaporfra)
                flag=1
                break
        """   
        if(flag==1):
            for i in np.arange(80, 200, 0.1):
                mix.T = float(i)
                mix.equalconstant = temp
                mix.solve(False)
                if(mix.vaporfra < 0.999999 and mix.vaporfra>2e-7):
                    print(i)
                    print(mix.vaporfra)
        """        
        if(flag==1):
            #print(ip)
            if(ip>tp):
                tp=ip
        else:
            break
    print(tp,ix)
    p.append(tp)
    x.append(ix)
plt.plot(x,p)
plt.savefig("./src/thermophysicalModels/pythoncode/plt.png")
print(p)
headers=["CH4","Pc"]
a=[x,p]
b=[[row[i] for row in a] for i in range(len(a[0]))]
with open('./src/thermophysicalModels/pythoncode/test.csv','w')as f:
    f_csv = csv.writer(f)
    f_csv.writerow(headers)
    f_csv.writerows(b)