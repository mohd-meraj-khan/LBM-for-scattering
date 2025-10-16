import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import scipy.special as sc




directory = 'data'
if not os.path.exists(directory):
    os.makedirs(directory)




ratio = np.arange(0.02, 50 + 0.0001, 0.001)



##n = 50




Fx_avg = np.zeros(len(ratio))



for i in range(len(ratio)):

    if (ratio[i] < 0.1):
        n = 20
    elif (ratio[i] < 1):
        n = 50
    elif (ratio[i] < 10):
        n = 100
    elif (ratio[i] < 20):
        n = 200
    elif (ratio[i] < 50):
        n = 350

    k1a = 2 * np.pi * ratio[i]
    
    sum = 0
    
    for l in range(-n, n+1, 1):

        fx = np.real(((0-1j)**l / sc.hankel2(l, k1a)) * np.conjugate(((0-1j)**(l+1) / sc.hankel2(l+1, k1a))))

        sum = sum + fx


    Fx_avg[i] = - 1 / (2 * np.pi**3 * ratio[i]) * np.real(sum)



np.save(directory+"/Fx_PEC_Exact.txt", Fx_avg)



########################################################################################################################################################
########################################################################################################################################################





ratio = [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3,
                0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50]





Fx_avg = np.zeros(len(ratio))



for i in range(len(ratio)):

    if (ratio[i] < 0.1):
        n = 20
    elif (ratio[i] < 1):
        n = 50
    elif (ratio[i] < 10):
        n = 100
    elif (ratio[i] < 20):
        n = 200
    elif (ratio[i] < 50):
        n = 350

    k1a = 2 * np.pi * ratio[i]
    
    sum = 0
    
    for l in range(-n, n+1, 1):

        fx = np.real(((0-1j)**l / sc.hankel2(l, k1a)) * np.conjugate(((0-1j)**(l+1) / sc.hankel2(l+1, k1a))))

        sum = sum + fx


    Fx_avg[i] = - 1 / (2 * np.pi**3 * ratio[i]) * np.real(sum)



np.save(directory+"/Fx_PEC_Exact_Err.txt", Fx_avg)



########################################################################################################################################################
########################################################################################################################################################
		
sys.exit()
###################################


plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)

################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.05), dpi=600, constrained_layout = True)




plt.plot(ratio, Fx_avg, 'k-')


plt.xlabel(r'$\phi (^o)$')
plt.ylabel(r'$\sigma / \lambda$')

plt.grid(which='both')


plt.xscale('log')

    
plt.savefig(directory+'/Fx_PEC1.svg')
plt.close()
################################################################################
