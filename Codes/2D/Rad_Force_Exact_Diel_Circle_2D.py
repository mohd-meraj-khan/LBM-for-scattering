import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import scipy.special as sc





directory = 'data'
if not os.path.exists(directory):
    os.makedirs(directory)


ratio = np.arange(0.02, 0.1+0.000001, 0.00001)



n = 70

er = 2


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


np.save(directory+"/Fx_Exact_er_{}.npy".format(er), Fx_avg)






###################################

ratio = [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]

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


np.save(directory+"/Fx_Exact_er_{}_Err.npy".format(er), Fx_avg)






###################################


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


##plt.yscale('log')

    
plt.savefig(directory+'/Fx_er_{}.svg'.format(er))
plt.close()
################################################################################

