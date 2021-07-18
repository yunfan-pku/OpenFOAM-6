import VLE
from VLE import DoubleVector as vec
import numpy as np
import matplotlib.pyplot as plt
import csv
import matplotlib.pyplot as plt
R = 1.38065e-23*6.0221417930e+23
mix = VLE.solver_new(
    "/home/zhy/Documents/git/OpenFOAM-6/src/thermophysicalModels/newtest/")
mix.m_specie = ["O2", "CH4"]
mix.reset()
mix.comp = [0.5, 0.5]
mix.P = 1e5
#mix.P = 1e5
temp = []
vtemp = []
ex = []


for T in np.arange(90, 110, 0.1):
    mix.T = float(T)
    mix.TPn_flash()
    temp.append(T)
    vtemp.append(mix.vaporfra)
    ex.append(mix.comp_gas[0])

    #print(T, mix.vaporfra)
plt.figure()
plt.plot(temp, vtemp)
font1 = {'family': 'DejaVu Sans', 'weight': 'normal', 'size': 14, }
plt.xlabel("T(K)", font1)
plt.ylabel("Gas mole fraction", font1)
plt.savefig("TPD2.png", dpi=300)
plt.close()

plt.figure()
plt.plot(temp, ex)
font1 = {'family': 'DejaVu Sans', 'weight': 'normal', 'size': 14, }
plt.xlabel("T(K)", font1)
plt.ylabel("O2 mole fraction_gas", font1)
plt.savefig("gasO2.png", dpi=300)
plt.close()


mix.T = 100
x = []
G = []
dGl = []
mix.comp = [0.0, 1.0]
Gl = mix.G_Mole()
print("Gl", Gl)
mix.comp = [1.0, 0.0]
Gr = mix.G_Mole()
TPD = []
print("Gr", Gr)
for x1 in np.arange(0, 1.01, 0.001):
    x.append(float(x1))
    mix.comp = [float(x1), 1-float(x1)]
    # G.append(mix.Gibbs_single()*mix.W()-float(x1)*Gr-(1-float(x1))*Gl)
    G.append(mix.G_Mole()-float(x1)*Gr-(1-float(x1))*Gl)
    # G.append(mix.W())
    mix.Ln_fugacityCoefficient()
    comp = mix.comp
    lnphi = mix.ret
    #dG = R*mix.T*(np.log(comp[0])+lnphi[0]-np.log(comp[1])-lnphi[1])
    dGl.append((R*mix.T*lnphi[0]+mix.muideal_Mole(0)) -
               (R*mix.T*lnphi[1]+mix.muideal_Mole(1))-(Gr-Gl))
    if mix.solveTPD_BFGS():
        TPD.append(0)
    else:
        TPD.append(-100)

mix.T = 100
mix.P = 1e5
x1 = 0.24
mix.comp = [x1, 1-x1]
mix.Ln_fugacityCoefficient()
lnphi = mix.ret
dGdx = (R*mix.T*lnphi[0]+mix.muideal_Mole(0)) - \
    (R*mix.T*lnphi[1]+mix.muideal_Mole(1))-(Gr-Gl)
G0 = mix.G_Mole()-float(x1)*Gr-(1-float(x1))*Gl

xt = np.linspace(0, 1.0, 21)
yt = dGdx*(xt-x1)+G0
print("xt=", xt)
print("yt=", yt)
print("dGdx=", dGdx)
print("TPD", TPD)
mix.TPn_flash()

print("gas", mix.comp_gas)
print("liq", mix.comp_liq)
print("vf", mix.vaporfra)

plt.plot(x, G)
plt.plot(xt, yt)
plt.plot(x, TPD)
plt.vlines(0.7280245714072449, -250, 0, colors="r", linestyles="dashed")
plt.vlines(0.23988303926952834, -250, 0, colors="r", linestyles="dashed")
#plt.hlines(234, 0, 1, colors="r", linestyles="dashed")
plt.savefig("G2.png", dpi=300)
plt.close()
plt.plot(x, dGl)
plt.savefig("dG2.png", dpi=300)
mix.T = 243


mix.T = 100
mix.P = 1e5
mix.comp = [0.2, 0.8]


mix.Ln_fugacityCoefficient()
comp = mix.comp
lnphi = mix.ret
# dG=R*mix.T*(np.log(comp[0])+lnphi[0]-np.log(comp[1])-lnphi[1])
dG = R*mix.T*lnphi[0]+mix.muideal_Mole(0)
# dG=mix.muideal_Mole(0)
Gideal1 = mix.Gideal()*mix.W()*(1)

a1 = mix.A_single()*mix.W()*1
g1 = mix.G_Mole()
G1 = mix.Gibbs_single()*mix.W()
dgdx = (R*mix.T*lnphi[0]+mix.muideal_Mole(0)-G1)/(2*0.8)
v1 = mix.z_single()*1.0*R*mix.T/mix.P
dx = 0.0001
mix.comp = [(0.2+dx)/(1+dx), 0.8/(1+dx)]
g2 = mix.G_Mole()*(1+dx)
G1 = mix.Gibbs_single()*mix.W()
Gideal2 = mix.Gideal()*mix.W()*(1+dx)
print((g2-g1)/dx, dG)
mix.P = mix.P+43075.8072857
v2 = mix.z_single()*(1.0+dx)*R*mix.T/mix.P
a2 = mix.A_single()*mix.W()*(1+dx)

print(v1, v2, v2-v1)


print((a2-a1)/dx, (g2-g1)/dx, dG-(g2-g1)/dx)

mix.T = 100
mix.P = 1e5
x1 = 0.6
mix.comp = [x1, 1-x1]
mix.Ln_fugacityCoefficient()
lnphi = mix.ret
G1 = mix.G_Mole()
Gi1 = mix.Gideal_Mole()
dgdx = (R*mix.T*lnphi[0]+mix.muideal_Mole(0)-G1)/(0.8)
dgdx = (R*mix.T*lnphi[0]+mix.muideal_Mole(0)) - \
    (R*mix.T*lnphi[1]+mix.muideal_Mole(1))

dgidealdx = mix.muideal_Mole(0)-mix.muideal_Mole(1)
dx = 1e-6

dgdx1 = R*mix.T*lnphi[0]+mix.muideal_Mole(0)
mix.comp = [(x1+dx)/(1+dx), (1-x1)/(1+dx)]
G2 = mix.G_Mole()*(1+dx)
Gi2 = mix.Gideal_Mole()
print((G2-G1)/dx, dgdx1, (G2-G1)/dx-dgdx1)

dx2 = -dx

dgdx2 = R*mix.T*lnphi[1]+mix.muideal_Mole(1)
mix.comp = [(x1+dx)/(1+dx+dx2), (1-x1+dx2)/(1+dx+dx2)]
G3 = mix.G_Mole()*(1+dx+dx2)

dgdx2 = R*mix.T*lnphi[1]+mix.muideal_Mole(1)
print((G3-G2)/dx2, dgdx2, (G3-G2)/dx2-dgdx2)


print((G3-G1)/dx, dgdx1 + dgdx2*dx2/dx, dgdx,
      (G3-G1)/dx-(dgdx1 + dgdx2*dx2/dx))


mix.T = 100
mix.P = 1e5
x1 = 0.2
mix.comp = [x1, 1-x1]
mix.Ln_fugacityCoefficient()
lnphi = mix.ret
G1 = mix.G_Mole()-float(x1)*Gr-(1-float(x1))*Gl
dgdx = (R*mix.T*lnphi[0]+mix.muideal_Mole(0)) - \
    (R*mix.T*lnphi[1]+mix.muideal_Mole(1))-(Gr-Gl)

dx = 1e-8
mix.comp = [(x1+dx), (1-x1-dx)]
G2 = mix.G_Mole()-float(x1+dx)*Gr-(1-float(x1+dx))*Gl


print((G2-G1)/dx, dgdx, (G2-G1)/dx-dgdx)
print(G1)
