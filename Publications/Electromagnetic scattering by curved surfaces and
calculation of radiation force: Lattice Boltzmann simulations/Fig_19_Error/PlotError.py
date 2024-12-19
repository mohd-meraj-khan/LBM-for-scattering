import numpy as np
import matplotlib.pyplot as plt

from scipy.signal import find_peaks
from scipy.optimize import curve_fit




fx_91 = np.array([0.5427, 0.5414, 0.5301, 0.5253, 0.5248, 0.5239, 0.5255, 0.5215, 0.5213, 0.5203, 0.5248, 0.5245, 0.5266, 0.5246, 0.5243, 0.5254])


fx_91 = np.loadtxt('Fx_Avg_ratio_0.91.txt')

a_91 = np.array([25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175])

fx_exact_91 = 0.5219

fx_err_91 = (fx_91 / fx_exact_91 - 1)*100



fx_93 = np.array([0.6361, 0.6118, 0.6161, 0.6187, 0.6348, 0.6288, 0.6248, 0.6248, 0.6229, 0.6225, 0.6290, 0.6266, 0.6249, 0.6231, 0.6260, 0.6255])


fx_93 = np.loadtxt('Fx_Avg_ratio_0.93.txt')

a_93 = np.array([25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175])

fx_exact_93 = 0.6209

fx_err_93 = (fx_93 / fx_exact_93 - 1)*100




fx_95 = np.array([1.4242, 0.9306, 1.0670, 1.0066, 1.4508, 1.3285, 1.1766, 1.1760, 1.1946, 1.1941, 1.4446, 1.3360, 1.3649,
               1.3544, 1.4169, 1.3977, 1.4280, 1.4022, 1.4103, 1.4017, 1.3822, 1.4132, 1.4226, 1.4037, 1.4318, 1.3979, 1.4182, 1.4119, 1.4150, 1.4178, 1.4096])

fx_95 = np.loadtxt('Fx_Avg_ratio_0.95.txt')

a_95 = np.array([25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175])

fx_exact_95 = 1.4465

fx_err_95 = (fx_95 / fx_exact_95 - 1)*100





############################################################################################################################################################################



plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)



##################################################################################################

fig, ax = plt.subplots(figsize = (3.35, 2), dpi=600, constrained_layout = True)

ax.plot(a_91, np.abs(fx_err_91), 'k-o')
ax.plot(a_93, np.abs(fx_err_93), 'r-o')
ax.plot(a_95, np.abs(fx_err_95), 'b-o')

ax.legend([r'$a/ \lambda = 0.91$', r'$a/ \lambda = 0.93$', r'$a/ \lambda = 0.95$'])

ax.set_ylabel(r'$| \frac{\left< F_{LBM} \right>}{\left< F_{exact} \right>} - 1 | \times 100$')
ax.set_xlabel(r'$a/ \Delta x$')

plt.savefig('Fx_Mie_err1.svg')
plt.close(fig)



##################################################################################################

print(find_peaks(fx_err_95))

X = np.array([a_95[1], a_95[3], a_95[7], a_95[9], a_95[11], a_95[13], a_95[15], a_95[17], a_95[20], a_95[23], a_95[25], a_95[27], a_95[29]])
Y = np.array([np.abs(fx_err_95[1]), np.abs(fx_err_95[3]), np.abs(fx_err_95[7]), np.abs(fx_err_95[9]), np.abs(fx_err_95[11]), np.abs(fx_err_95[13]), np.abs(fx_err_95[15]), np.abs(fx_err_95[17]), 
        np.abs(fx_err_95[20]), np.abs(fx_err_95[23]), np.abs(fx_err_95[25]), np.abs(fx_err_95[27]), np.abs(fx_err_95[29])])


def curve_func(x, a, b):
    return a * np.exp(-b * x)

initial_guess = [0.6, 0.19]
fit_params, _ = curve_fit(curve_func, X, Y, p0=initial_guess)
a, b = fit_params
print(a, b)


x = np.arange(30, 170, 1)
y = a * np.exp(-b * x)



##################################################################################################

fig, ax = plt.subplots(figsize = (3.35, 2), dpi=600, constrained_layout = True)

ax.plot(X, Y, 'ko')


ax.set_ylabel(r'$| \frac{\left< F_{LBM} \right>}{\left< F_{exact} \right>} - 1 | \times 100$')
ax.set_xlabel(r'$a/ \Delta x$')

plt.yscale('log')
plt.xscale('log')

plt.savefig('Fx_Mie_peaks_0.951.svg')
plt.close(fig)



##################################################################################################




