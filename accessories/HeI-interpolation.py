#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  9 15:59:03 2025

@author: hmonteiro
"""
import numpy as np
import matplotlib.pyplot as plt
import sys, os


temp=np.linspace(40,20000,100)
T4 = temp/1e4


plt.figure()
plt.xscale('log')
plt.yscale('log')
plt.ylim(0,100)
#plt.ylim(temp.min(),temp.max())

denint = 1
NeUsed = 8000.
x1=1.73*(T4**(-1.091))*np.exp(-0.024/T4)
x2=0.842*(T4**(0.304))*np.exp(0.821/T4)

HeIRecLines = x1+((x2-x1)*(NeUsed-100.**denint)/(100.**(denint+1)-100.**(denint)))

plt.plot(temp,HeIRecLines, 'C0')
plt.plot(temp,x1, '.C0',alpha=0.2)
plt.plot(temp,x2, '.C0',alpha=0.2)


denint = 1
NeUsed = 1000.
x1=1.73*(T4**(-1.091))*np.exp(-0.024/T4)
x2=0.842*(T4**(0.304))*np.exp(0.821/T4)

HeIRecLines = x1+((x2-x1)*(NeUsed-100.**denint)/(100.**(denint+1)-100.**(denint)))

plt.plot(temp,HeIRecLines, 'C1')
plt.plot(temp,x1, '.C1',alpha=0.2)
plt.plot(temp,x2, '--C1',alpha=0.2)


# power law extrapolation

# slope from Te=5000 and 6000
x1_5000 = 1.73*((0.5)**(-1.091))*np.exp(-0.024/0.5)
x1_6000 = 1.73*((0.6)**(-1.091))*np.exp(-0.024/0.6)

x2_5000 = 0.842*((0.5)**(0.304))*np.exp(0.821/0.5)
x2_6000 = 0.842*((0.6)**(0.304))*np.exp(0.821/0.6)

m1 = (np.log10(x1_5000)-np.log10(x1_6000))/(np.log10(0.5)-np.log10(0.6))
m2 = (np.log10(x2_5000)-np.log10(x2_6000))/(np.log10(0.5)-np.log10(0.6))

plt.plot(temp,(x1_5000/0.5**m1)*T4**m1,'C2')
plt.plot(temp,(x2_5000/0.5**m2)*T4**m2,'C3')





