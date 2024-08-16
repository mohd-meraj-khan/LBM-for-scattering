import sys



# number of parallel threads
N = 1


# Accessing command line arguments
parameters = sys.argv


######################


a, ratio = 25, 0.5   # ratio = a / wavelength
n = 4
Nx, Ny = n*a, n*a  # size of the computational domain

er1, mur1, er2, er3 = 1, 1, 2, 5   # material properties i.e. permittivity and permeabilty


#############################
# boundary of EM wave source
xloc = 0
ymin = 0
ymax = Ny
#############################


noOfPeriods = 5
noOfReflections = 0



