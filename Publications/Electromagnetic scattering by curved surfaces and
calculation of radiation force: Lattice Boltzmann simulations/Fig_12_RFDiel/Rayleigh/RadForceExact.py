import numpy as np ; from numpy.linalg import *
import sys
import os
import csv
import math
import scipy.special as sc
from math import e
import cmath
import matplotlib.pyplot as plt




ratio = np.arange(0.02, 0.1+0.0001, 0.001)

n = 70

er = 5


Fx_avg = np.zeros(len(ratio))

print(len(ratio))


for i in range(len(ratio)):


    k1a = 2 * np.pi * ratio[i]
    k2a = 2 * np.pi * ratio[i] * np.sqrt(er)

    sum = 0
    
    for l in range(0, n+1, 1):

        alphan  = np.sqrt(er) * sc.hankel2(l, k1a) * sc.jvp(l, k2a) - sc.h2vp(l, k1a) * sc.jv(l, k2a)
        alphan1 = np.sqrt(er) * sc.hankel2(l+1, k1a) * sc.jvp(l+1, k2a) - sc.h2vp(l+1, k1a) * sc.jv(l+1, k2a)

        fx = 1 / (np.pi**5 * ratio[i]**2) * (er - 1)**2 * (sc.jv(l, k2a)**2 * sc.jv(l+1, k2a)**2) / np.abs(alphan * alphan1)**2

##        if (fx <= 1e-10):
####            print(l)
##            break

        sum = sum + fx


    Fx_avg[i] = sum


Fx = open("FxExact_er_{}.txt".format(er), "w")
np.savetxt(Fx, Fx_avg)
Fx.close()

