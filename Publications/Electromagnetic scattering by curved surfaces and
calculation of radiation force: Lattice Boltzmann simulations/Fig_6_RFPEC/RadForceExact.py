import numpy as np

import sys
import os
import scipy.special as sc
from math import e
import cmath



directory = 'data'
if not os.path.exists(directory):
    os.makedirs(directory)


ratio = np.arange(0.02, 0.1 + 0.0001, 0.001)

##ratio = np.arange(0.5, 4.0 + 0.0001, 0.01)

n = 50




Fx_avg = np.zeros(len(ratio))



for i in range(len(ratio)):

    k1a = 2 * np.pi * ratio[i]
    
    sum = 0
    
    for l in range(-n, n+1, 1):

        fx = np.real(((0-1j)**l / sc.hankel2(l, k1a)) * np.conjugate(((0-1j)**(l+1) / sc.hankel2(l+1, k1a))))

        sum = sum + fx


    Fx_avg[i] = - 1 / (2 * np.pi**3 * ratio[i]) * np.real(sum)



Fx = open(directory+"/FxRayleighExact.txt", "w")

##Fx = open(directory+"/FxMieExact.txt", "w")

np.savetxt(Fx, Fx_avg)
Fx.close()


########################################################################################################################################################
########################################################################################################################################################
			

