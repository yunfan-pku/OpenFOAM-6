import VLE
from VLE import DoubleVector as vec
import numpy as np
import matplotlib.pyplot as plt
import csv
import matplotlib.pyplot as plt
mix = VLE.solver_new(
    "/home/zhy/Documents/git/OpenFOAM-6/src/thermophysicalModels/newtest/")
mix.m_specie = ["CO2", "H2O", "CH4"]
mix.reset()
mix.comp = [0.199932, 0.000339, 0.799729]

"""
mix.P = 0.15e7
mix.T = 175
ff =[]
for p in np.arange(1e5,20e5,1e5):
    mix.P =p
    mix.TPn_flash_New_TPD()
    ff.append(str(p)+":"+str(mix.vaporfra))
print(ff)





"""
#mix.P = 1e5
temp = []
vtemp = []
TPD = []
Tlist = np.linspace(100, 300, 101)
Plist = np.linspace(1e5, 100e5, 100)
T, P = np.meshgrid(Tlist, Plist)
TPDout = T*T
vf = T*T



for t in range(len(Tlist)):
    for p in range(len(Plist)):
        mix.T = float(T[p, t])
        mix.P = float(P[p, t])
        #print("start!",mix.T,mix.P)
        mix.TPn_flash_New_TPD()
        #print("end!",mix.T,mix.P)
        vf[p,t]=mix.vaporfra
        # temp.append(T)
        # vtemp.append(mix.vaporfra)
        if mix.solveTPD_BFGS():
            TPDout[p, t] = 1.0
        else:
            TPDout[p, t] = 0.0
        print(mix.T,mix.P, mix.solveTPD_BFGS())
        

plt.contourf(T, P, TPDout)

plt.savefig("TPDout.png", dpi=300)
plt.close()
#plt.contourf(T, P, vf,levels=[0.0001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99999])
plt.contourf(T, P, vf,levels=[0.00001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.96,0.97,0.98,0.99,0.99999])
plt.colorbar();
#plt.scatter([175,175],[21e5,22e5],color='r',s=5)
plt.savefig("TPDout2_new2.png", dpi=300)



"""
plt.plot(temp, vtemp)
plt.plot(temp, TPD)
font1 = {'family': 'Times New Roman', 'weight': 'normal', 'size': 14, }
plt.xlabel("T(K)", font1)
plt.ylabel("Gas mole fraction", font1)
plt.savefig("TPD1_org_0.png", dpi=300)
mix.T = 244
mix.TPn_flash()
print(mix.comp_liq, mix.comp_gas, mix.vaporfra)
mix.fugacityCoefficient(1, vec(mix.comp_gas))
print(mix.ret)
mix.fugacityCoefficient(0, vec(mix.comp_liq))
print(mix.ret)
print(mix.solveTPD_BFGS())
"""
