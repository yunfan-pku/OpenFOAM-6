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
mix.comp = [0.425, 0.075, 0.5]
mix.P = 30000000
#mix.P = 1e5
temp = []
vtemp = []
TPD = []
for T in np.arange(20, 600, 1):
    mix.T = float(T)
    #mix.TPn_flash()
    mix.TPn_flash_New()
    temp.append(T)
    vtemp.append(mix.vaporfra)
    if mix.solveTPD_BFGS():
        TPD.append(0.5)
    else:
        TPD.append(0.3)

    print(T, mix.vaporfra, mix.solveTPD_BFGS())
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
