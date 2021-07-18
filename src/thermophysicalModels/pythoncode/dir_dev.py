import VLE
from VLE import DoubleVector as vec
import numpy as np
import matplotlib.pyplot as plt
import csv
mix = VLE.solver_new(
    "/home/zhy/Documents/git/OpenFOAM-6/src/thermophysicalModels/newtest/")
"""
mix.T = 270
mix.P = 30000000
mix.comp = [0.5, 0.5]

mix.TPn_flash()
vf1 = mix.vaporfra
# print(vf1)
# print(mix.comp_liq)
# print(mix.comp_gas)
# mix.comp = mix.comp_gas
comp_gas = mix.comp_gas

comp_gas = [0.4591236483139583, 0.0007217431022268374, 0.5401546085838149]
comp_liq = [9.892031327512852e-05, 0.9998985121031675, 2.5675835573989294e-06]
# mix.comp = comp_gas


dzdt = mix.dZdT(1)
dadt = mix.dAdT()
dbdt = mix.dBdT()
cp = mix.Cp()
ha1 = mix.Ha()
mix.ddT_Ln_fugacityCoefficient(1)
dfcdT_gas = np.array(mix.ret)

# mix.comp = comp_liq
mix.dvidT()
dbdP = np.sum(np.array(mix.ret))
print(dbdP)

mix.ddT_Ln_fugacityCoefficient(0)
dfcdT_liq = np.array(mix.ret)

z1 = mix.Z(1)
a1 = mix.A()
b1 = mix.B()
mix.dvidP()
mix.Ln_fugacityCoefficient(1)
ln_gas_g_f1 = np.array(mix.ret)
# print(mix.ret)
# mix.comp = comp_liq
mix.Ln_fugacityCoefficient(1)
ln_liq_g_f1 = np.array(mix.ret)


dT = 1e-5
dP = 1

mix.T = mix.T+dT
# mix.P = mix.P+dP
mix.TPn_flash()
vf2 = mix.vaporfra
ha2 = mix.Ha()
# print(vf2)
"""
"""
mix.comp = comp_gas
# print(z1*z1*z1+(b1-1)*z1*z1+(a1-2*b1-3*b1*b1)*z1+(b1*b1*b1+b1*b1-a1*b1))
dT = 1e-6
dxi = 1e-6
comp_gas[0] += dxi
# mix.T = mix.T+dT
mix.comp = comp_gas
z2 = mix.Z(1)
a2 = mix.A()
b2 = mix.B()
mix.Ln_fugacityCoefficient(1)
ln_g_f2 = np.array(mix.ret)
# print((z2-z1)/dT)
# print((a2-a1)/dT)
# print((b2-b1)/dT)
# print(z2*z2*z2+(b2-1)*z2*z2+(a2-2*b2-3*b2*b2)*z2+(b2*b2*b2+b2*b2-a2*b2))
dzdtt = (z2-z1)/dT

# print(3*z1*z1*dzdt+(dbdt)*z1*z1+2*(b1-1)*z1*dzdt+(dadt-2*dbdt-3*2*b1*dbdt)
#      * z1+(a1-2*b1-3*b1*b1)*dzdt+(3*b1*b1*dbdt+2*b1*dbdt-dadt*b1-a1*dbdt))
# print(3*z1*z1*dzdtt+(dbdt)*z1*z1+2*(b1-1)*z1*dzdtt+(dadt-2*dbdt-3*2*b1*dbdt)
#      * z1+(a1-2*b1-3*b1*b1)*dzdtt+(3*b1*b1*dbdt+2*b1*dbdt-dadt*b1-a1*dbdt))
print((ln_g_f2-ln_gas_g_f1)/dT)
print(dfcdT)
print((ln_g_f2-ln_gas_g_f1)/dT-dfcdT)
# print(comp_gas*np.exp(ln_gas_g_f1))
# print(comp_liq*np.exp(ln_liq_g_f1))
"""
comp = [0.425, 0.075, 0.5]
comp_gas = [0.4591236483139583, 0.0007217431022268374, 0.5401546085838149]
comp_liq = [9.892031327512852e-05, 0.9998985121031675, 2.5675835573989294e-06]

mix.T = 270
mix.P = 30000000
mix.comp = comp


mix.TPn_flash()
dtdx0 = mix.dTdXi(0)
dtdp = mix.dTdP()
dtdha = mix.dTdH()
Z1 = mix.Z()
rho1 = mix.rho()
dZdt = mix.dZdT()
drhodt = mix.drhodT()
drhodp = mix.drhodP()
drhodx0 = mix.drhodXi(0)
drhodx0_HP = mix.drhodXi_HP(0)
dZdx0 = mix.dZdXi(0)
drhoPdP = mix.drhoPdP_HP()
drhoPdxi = mix.drhoPdXi_HP(0)
drhodxi = mix.drhodXi_HP(0)
drhodP = mix.drhodP_HP()
dvfdH = mix.dvfdH_HP()
dvfdP = mix.dvfdP_HP()
dvfdxi = mix.dvfdXi_HP(0)
print("Z1=", Z1)


N1_gas = mix.comp_gas
vf1 = mix.vaporfra
Ntotal1 = np.sum(np.array(comp))
print(Ntotal1)
v1 = vf1*np.array(mix.comp_gas)*Ntotal1
ha1 = mix.Ha()
ha1_gas = mix.Ha_singlePhase(0, vec(mix.comp_liq))
ygas1 = vf1*mix.W(vec(mix.comp_gas))/mix.W(vec(mix.comp))
w1 = mix.W()
phi1 = mix.P/mix.rho()

comp_liq1 = np.array(mix.comp_liq)
comp_gas1 = np.array(mix.comp_gas)
dhdx = mix.dHadXi(0)
dhdP = mix.dHadP()
dhdT = mix.Cp()
mix.dvidXi(0)
dvid0 = np.array(mix.ret)
mix.dvidXi(1)
dvid1 = np.array(mix.ret)
mix.dvidXi(2)
dvid3 = np.array(mix.ret)


dxi = 1e-7
dP = 1e-2

dT = dtdx0*dxi
print("dT=", dT)
comp2 = np.array(comp)
comp2[0] += dxi
# comp2[1] -= 3*dxi
# comp2[2] += 2*dxi
Ntotal2 = np.sum(comp2)

print(Ntotal2)
comp2 /= Ntotal2
#comp[1] -= dxi
print(comp2)
mix.comp = comp2
mix.T = mix.T+dT
#mix.P = mix.P+dP
mix.TPn_flash()
Z2 = mix.Z()
phi2 = mix.P/mix.rho()
rho2 = mix.rho()
ha2 = mix.Ha()
vf2 = mix.vaporfra
w2 = mix.W()
v2 = vf2*np.array(mix.comp_gas)*Ntotal2
ygas2 = vf2*mix.W(vec(mix.comp_gas))/mix.W(vec(mix.comp))
ha2_gas = mix.Ha_singlePhase(0, vec(mix.comp_liq))
comp_liq2 = np.array(mix.comp_liq)
comp_gas2 = np.array(mix.comp_gas)
# print((vf2-vf1)/dxi)

print(dvid0-dvid1)

print(np.sum(dvid0-dvid1))
# print(np.sum(dvid0-dvid1))
# print(dvid0)
# print(dvid1)


print((ha2_gas-ha1_gas)/dxi)
print((comp_liq2-comp_liq1)/dxi)
print((comp_gas2-comp_gas1)/dxi)
print((v2-v1)/dxi)
print(dvid0)

print(dT/(ha2-ha1))
print(dtdha)
print(dhdT)
print(dhdx)

print((Z2-Z1)/dxi)
print(dZdx0)
print((ha2-ha1)/dxi)
print((vf2-vf1)/(dxi))
print(dvfdxi)
