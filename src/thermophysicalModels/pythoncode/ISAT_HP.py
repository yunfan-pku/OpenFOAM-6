#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 21:53:30 2021

@author: zhy
"""


import VLE
from VLE import DoubleVector as vec
import numpy as np
import matplotlib.pyplot as plt
import csv
mix = VLE.solver_new(
    "/home/zhy/Documents/git/OpenFOAM-6/src/thermophysicalModels/newtest/")
mix2 = VLE.solver_new(
    "/home/zhy/Documents/git/OpenFOAM-6/src/thermophysicalModels/newtest/")
mix3 = VLE.solver_new(
    "/home/zhy/Documents/git/OpenFOAM-6/src/thermophysicalModels/newtest/")
mix.m_specie = ["CO2", "H2O"]
mix.reset()
mix.comp = [0.69995272395133235, 0.30004727604866765]
T = mix.T_HsP(226106.08046651725, 10079795.339350015, 550)
mix.T = T
mix.P = 10079795.339350015
mix.TPn_flash()
print(T)
print(mix.Hs())
v1 = np.array([0.69995272395133235, 0.30004727604866765,
               226106.08046651725, 10079795.339350015])


mix2.comp = [0.69995272395126051, 0.30004727604873943]
mix2.T = 577.69502777036882
mix2.T_HsP(220089.01111327758, 23000000, 550)
mix2.P = 12338482.855487322
mix2.TPn_flash()
phi = mix2.rho()/mix2.P
print(phi)
gr = [mix2.drhoPdXi_HsP(0), mix2.drhoPdXi_HsP(
    1), mix2.drhoPdH_HsP(), mix2.drhoPdP_HsP()]

v2 = np.array([0.69995272395133512, 0.30004727604866477,
               220089.01111327758, 10017243.890373562])

comp2 = np.array([0.69995272395133512, 0.30004727604866477])
dx = 1e-5
comp2[0] += dx
comp2 /= np.sum(comp2)
dH = (v1-v2)[2]
dP = 1e-1
mix3.comp = comp2
mix3.T = 577.69502777036882
mix3.T_HsP(220089.01111327758, 10017243.890373562, 550)
mix2.T_HsP(220089.01111327758, 10017243.890373562, 550)
mix3.P = 12338482.855487322
mix3.TPn_flash()
phib = mix3.rho()/mix3.P
print((phib-phi)/dx)
print((mix3.T-mix2.T)/dx)
print(mix2.dTdXi_HsP(0))
print((v1-v2)*gr)
print(mix3.vaporfra, mix2.vaporfra, (mix3.vaporfra-mix2.vaporfra)/dx)
print(mix2.dvfdXi_HsP(0))



mix2.drhoPdXHP_HsP()
print(gr)
print(mix2.ret)
print((np.array(mix2.ret)-np.array(gr))/gr)
print((mix3.vaporfra*mix3.W(vec(mix3.comp_gas))/mix3.W() -
       mix2.vaporfra*mix2.W(vec(mix2.comp_gas))/mix2.W())/dx)
