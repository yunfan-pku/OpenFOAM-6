import VLE
from VLE import DoubleVector as vec
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import csv
import matplotlib.pyplot as plt
mix = VLE.solver_new(
    "/home/zhy/Documents/git/OpenFOAM-6/src/thermophysicalModels/newtest/")
mix.m_specie = ["CO2", "H2O", "CH4"]
mix.reset()
R = 1.38065e-23*6.0221417930e+23
mix.comp = [0.199932, 0.000339, 0.799729]
mix.P = 800000
mix.T = 175
#mix.P = 1e5
temp = []
vtemp = []
TPD = []
TPDr = mix.solveTPD_BFGS()
print(TPDr)
x = []
y = []
G = []
mix.comp = [1.0, 0, 0]
G1 = mix.G_Mole()
mix.comp = [0.0, 1.0, 0.0]
G2 = mix.G_Mole()
mix.comp = [0.0, 0.0, 1.0]
G3 = mix.G_Mole()
for x1 in np.arange(0, 1, 0.01):
    for x2 in np.arange(0, 1-x1, 0.01):
        x.append(x1)
        y.append(x2)

        mix.comp = [float(x1), float(x2), float(1-x1-x2)]
    # G.append(mix.Gibbs_single()*mix.W()-float(x1)*Gr-(1-float(x1))*Gl)
        G.append(mix.G_Mole()-float(x1)*G1-(float(x2))*G2-(1-float(x1+x2))*G3)
    # G.append(mix.W())
mix.comp = [0.199932, 0.000339, 0.799729]
comp =[0.199932, 0.000339, 0.799729]
mix.TPn_flash_New()
print(mix.comp_liq,mix.comp_gas,mix.vaporfra)
comp_liq=mix.comp_liq
comp_gas=mix.comp_gas
mix.comp=comp_liq
G_l=mix.G_Mole()-comp_liq[0]*G1-comp_liq[1]*G2-comp_liq[2]*G3
mix.comp=comp_gas
G_g=mix.G_Mole()-comp_gas[0]*G1-comp_gas[1]*G2-comp_gas[2]*G3
print(G_l,G_g,G_g*mix.vaporfra+G_l*(1-mix.vaporfra))
mix.comp=comp
mix.TPn_flash()
print(mix.comp_liq,mix.comp_gas,mix.vaporfra)
comp_liq=mix.comp_liq
comp_gas=mix.comp_gas
mix.comp=comp_liq
G_l=mix.G_Mole()-comp_liq[0]*G1-comp_liq[1]*G2-comp_liq[2]*G3
mix.comp=comp_gas
G_g=mix.G_Mole()-comp_gas[0]*G1-comp_gas[1]*G2-comp_gas[2]*G3
print(G_l,G_g)

comp=[0.199932, 0.000339, 0.799729]
mix.Ln_fugacityCoefficient()
lnphi = mix.ret
G0 = mix.G_Mole()-comp[0]*G1-comp[1]*G2-comp[2]*G3
dG0 = (R*mix.T*lnphi[0]+mix.muideal_Mole(0)) - \
    (R*mix.T*lnphi[2]+mix.muideal_Mole(2))-(G1-G3)
dG1 = (R*mix.T*lnphi[1]+mix.muideal_Mole(1)) - \
    (R*mix.T*lnphi[2]+mix.muideal_Mole(2))-(G2-G3)
print(G0,dG0,dG1)
with open('/home/zhy/Documents/git/OpenFOAM-6/src/thermophysicalModels/pythoncode/TPDtest'+str(int(mix.T))+'.csv', 'w')as f:
    f_csv = csv.writer(f)
    f_csv.writerow(x)
    f_csv.writerow(y)
    f_csv.writerow(G)
    f_csv.writerow([G0,dG0,dG1])



#plt.savefig("TPDtest.png", dpi=300)
"""
plt.plot(temp, vtemp)
plt.plot(temp, TPD)
font1 = {'family': 'Times New Roman', 'weight': 'normal', 'size': 14, }
plt.xlabel("T(K)", font1)
plt.ylabel("Gas mole fraction", font1)
plt.savefig("TPD1.png", dpi=300)
mix.T = 243
mix.TPn_flash()
print(mix.comp_liq, mix.comp_gas, mix.vaporfra)
mix.fugacityCoefficient(1, vec(mix.comp_gas))
print(mix.ret)
mix.fugacityCoefficient(0, vec(mix.comp_liq))
print(mix.ret)
"""
# %matplotlib qt
