import os
from ctypes import CDLL, POINTER
from ctypes import c_size_t, c_double, c_float, c_int
import numpy as np


###############################################################################################################
########                                           SHARED LIBRARY                                   ###########
###############################################################################################################

# loading the shared file (c library)
path = os.getcwd()
myclib = CDLL(os.path.join(path, "LBM_Sequential.so"))

# defining 2D and 3D pointers (LBM runs in C, for that pointer is needed)
P2D = np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags="C")
P3D = np.ctypeslib.ndpointer(dtype=np.float32, ndim=3, flags="C")

# calculation of macroscopic fields (FUNCTION PROTOTYPE)
myclib.macroField.argtypes = [P3D, P2D, P2D, c_int, c_int, c_int]
myclib.macroField.restype  = None

# initilization of macroscopic fields (FUNCTION PROTOTYPE)
myclib.initializeField.argtypes = [P2D, P2D, P2D, P2D, P2D, P2D, c_int, c_int]
myclib.initializeField.restype  = None

# collision and streaming (FUNCTION PROTOTYPE)
myclib.collStream.argtypes = [P3D, P3D, P3D, P3D, P3D, P3D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, c_int, c_int, c_int]
myclib.collStream.restype  = None

myclib.collStreamForcingNode.argtypes = [P3D, P3D, P3D, P3D, P3D, P3D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, c_int, c_int, c_int, c_int, c_int, c_int]
myclib.collStreamForcingNode.restype  = None

###############################################################################################################
