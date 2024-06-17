import numpy as np



#################################################################

def carToPolar(r, phi, Nz=1, Ny=1, Nx=1, cy=0, cx=0):
    for i in range(Ny):
        for j in range(Nx):
            r[0, i, j] = np.sqrt((i - cy)**2 + (j - cx)**2)
            phi[0, i, j] = np.arctan2((i - cy), (j - cx))   
#################################################################



#################################################################

def circle(r, phi, er, er2, a, Nz=1, Ny=1, Nx=1):
    for i in range(Ny):
        for j in range(Nx):
            if (r[0, i, j] <= a):
                er[0, i, j]  = er2
#################################################################




#################################################################

def square(er, er2, a, Nz=1, Ny=1, Nx=1, cy=0, cx=0):
    for i in range(Ny):
        for j in range(Nx):
            Y = i - cy
            X = j - cx
            if (X <= a and X >= -a and Y <= a and Y >= -a):
                er[0, i, j]  = er2
#################################################################



#################################################################

def hexagon(er, er2, a, Nz=1, Ny=1, Nx=1, cy=0, cx=0):
    for i in range(Ny):
        for j in range(Nx):
            Y = i - cy
            X = j - cx
            if (Y <= X + 2*a and Y >= - X - 2*a and Y >= X - 2*a and Y <= -X + 2*a and X >= - a and X <= a):
                er[0, i, j]  = er2
#################################################################   




#################################################################

def corrugatedEllipse(r, phi, er, er2, a, Nz=1, Ny=1, Nx=1, A=1, N=0, epsilon=0):
    for i in range(Ny):
        for j in range(Nx):
            if (r[0, i, j] <= a * (1 / ((np.cos(phi[0, i, j]))**2 + (A * np.sin(phi[0, i, j]))**2) + epsilon * np.cos(N * phi[0, i, j]))):
                er[0, i, j]  = er2
#################################################################

                


#################################################################

def JanusHalf(r, er, a, theta, er2, er3, Nz=1, Ny=1, Nx=1, cy=0, cx=0):
    for i in range(Ny):
        for j in range(Nx):
            Y = i - cy
            X = j - cx
            if (theta == 0 or theta == 360):
                if (r[0, i, j] <= a and X <= 0):
                    er[0, i, j]  = er2
                elif (r[0, i, j] <= a and X > 0):
                    er[0, i, j]  = er3
            elif (theta == 180):
                if (r[0, i, j] <= a and X >= 0):
                    er[0, i, j]  = er2
                elif (r[0, i, j] <= a and X < 0):
                    er[0, i, j]  = er3
            elif (0 < theta < 180):
                if (r[0, i, j] <= a and Y >= X * np.tan(np.pi/2 - theta * np.pi/180)):
                    er[0, i, j]  = er2
                elif (r[0, i, j] <= a and Y < X * np.tan(np.pi/2 - theta * np.pi/180)):
                    er[0, i, j]  = er3
            elif (180 < theta < 360):
                if (r[0, i, j] <= a and Y <= X * np.tan(3*np.pi/2 - theta * np.pi/180)):
                    er[0, i, j]  = er2
                elif (r[0, i, j] <= a and Y > X * np.tan(3*np.pi/2 - theta * np.pi/180)):
                    er[0, i, j]  = er3
#################################################################




#################################################################

def JanusSegment(r, er, a, er2, er3, Nz=1, Ny=1, Nx=1, cy=0, cx=0, shift=0):
    for i in range(Ny):
        for j in range(Nx):
            Y = i - cy
            X = j - cx
            if (r[0, i, j] <= a and X <= -(a-shift)):
                er[0, i, j]  = er2
            elif (r[0, i, j] <= a and X > -(a-shift)):
                er[0, i, j]  = er3
#################################################################




#################################################################

def JanusSector(r, er, a, er2, er3, Nz=1, Ny=1, Nx=1, cy=0, cx=0, theta=0):
    for i in range(Ny):
        for j in range(Nx):
            Y = i - cy
            X = j - cx
            if (r[0, i, j] <= a):
                if (theta <= 90):
                    if (X <= 0 and np.tan(theta * np.pi/180)*X <= Y <= -np.tan(theta * np.pi/180)*X):
                        er[0, i, j]  = er2
                    else :
                        er[0, i, j]  = er3
                else :
                    if (X > 0 and -np.tan((180-theta) * np.pi/180)*X <= Y <= np.tan((180-theta) * np.pi/180)*X):
                        er[0, i, j]  = er3
                    else :
                        er[0, i, j]  = er2
#################################################################














