import VLE
from VLE import DoubleVector as vec
import numpy as np
import matplotlib.pyplot as plt
import csv
import matplotlib.pyplot as plt
mix = VLE.solver_new(
    "/home/zhy/Documents/git/OpenFOAM-6/src/thermophysicalModels/newtest/")
mix.m_specie = ["CO2","H2O"]
mix.reset()
mix.T = 473
mix.comp = [0.199932, 0.000339, 0.799729]
dx=1e-2
l=[]
lp=[]


for p in np.arange(1e6,60e6,1e6):
    for x in np.arange(0,1.0-dx,dx):
        mix.P=p
        mix.comp=[float(x),float(1-x)]
        s1=mix.solveTPD_BFGS()
        mix.comp=[float(x+dx),float(1-x-dx)]
        s2=mix.solveTPD_BFGS()
        if(s1!=s2):
            l.append(x)
            lp.append(p)
mix.T = 523
for p in np.arange(1e6,60e6,1e6):
    for x in np.arange(0,1.0-dx,dx):
        mix.P=p
        mix.comp=[float(x),float(1-x)]
        s1=mix.solveTPD_BFGS()
        mix.comp=[float(x+dx),float(1-x-dx)]
        s2=mix.solveTPD_BFGS()
        if(s1!=s2):
            l.append(x)
            lp.append(p)
mix.T = 548
for p in np.arange(1e6,60e6,1e6):
    for x in np.arange(0,1.0-dx,dx):
        mix.P=p
        mix.comp=[float(x),float(1-x)]
        s1=mix.solveTPD_BFGS()
        mix.comp=[float(x+dx),float(1-x-dx)]
        s2=mix.solveTPD_BFGS()
        if(s1!=s2):
            l.append(x)
            lp.append(p)
mix.T = 573
for p in np.arange(1e6,60e6,1e6):
    for x in np.arange(0,1.0-dx,dx):
        mix.P=p
        mix.comp=[float(x),float(1-x)]
        s1=mix.solveTPD_BFGS()
        mix.comp=[float(x+dx),float(1-x-dx)]
        s2=mix.solveTPD_BFGS()
        if(s1!=s2):
            l.append(x)
            lp.append(p)
mix.T = 623
for p in np.arange(1e6,60e6,1e6):
    for x in np.arange(0,1.0-dx,dx):
        mix.P=p
        mix.comp=[float(x),float(1-x)]
        s1=mix.solveTPD_BFGS()
        mix.comp=[float(x+dx),float(1-x-dx)]
        s2=mix.solveTPD_BFGS()
        if(s1!=s2):
            l.append(x)
            lp.append(p)

plt.plot(l,lp,".")
plt.xlim(0,1)
plt.ylim(0,60e6)
plt.savefig("PX2.png", dpi=300)
