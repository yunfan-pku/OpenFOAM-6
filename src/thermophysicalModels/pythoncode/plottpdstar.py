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
mix.comp =  [0.425, 0.075, 0.5]
mix.P = 30000000
mix.T = 421
mix.Ln_fugacityCoefficient()
lnphi = np.array(mix.ret)
comp= np.array(mix.comp)
d=np.log(comp)+lnphi
print(d)




mix.comp = [0.005, 0.999, 0.005]
mix.P = 30000000
mix.T = 421
mix.Ln_fugacityCoefficient()
lnphi = np.array(mix.ret)
comp= np.array(mix.comp)
z=np.log(comp)+lnphi

tpd=comp.dot(z-d)
print(tpd)

lam=tpd+1
wt=np.exp(-tpd)
print(np.exp(-tpd))
print(1-wt+wt*np.log(wt)+wt*tpd)
