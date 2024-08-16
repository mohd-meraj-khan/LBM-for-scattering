import os
from ctypes import CDLL, POINTER
from ctypes import c_size_t, c_double, c_float, c_int
import numpy as np


###############################################################################################################
########                                           SHARED LIBRARY                                   ###########
###############################################################################################################

# loading the shared file (c library)
path = os.getcwd()
myclib = CDLL(os.path.join(path, "LBM.so"))

# defining 3D and 4D pointers (LBM runs in C, for that pointer is needed)
P2D = np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags="C")
P3D = np.ctypeslib.ndpointer(dtype=np.float32, ndim=3, flags="C")

# calculation of macroscopic fields (FUNCTION PROTOTYPE)
myclib.macroField.argtypes = [P3D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, c_int, c_int, c_int, c_int]
myclib.macroField.restype  = None


# collision + streaming (FUNCTION PROTOTYPE)
myclib.collNotForcingNode.argtypes = [P3D, P3D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, c_int, c_int, c_int, c_int]
myclib.collNotForcingNode.restype  = None


myclib.collForcingNode.argtypes = [P3D, P3D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, c_int, c_int, c_int, c_int, c_int, c_int, c_int]
myclib.collForcingNode.restype  = None


myclib.streaming.argtypes = [P3D, P3D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, P2D, c_int, c_int, c_int, c_int]
myclib.streaming.restype  = None

###############################################################################################################
