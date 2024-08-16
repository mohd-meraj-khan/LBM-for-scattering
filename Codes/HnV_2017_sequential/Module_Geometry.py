import numpy as np



#################################################################

def carToPolar(r, phi, Ny=10, Nx=10, cy=0, cx=0):
    for i in range(Ny):
        for j in range(Nx):
            r[i, j] = np.sqrt((i - cy)**2 + (j - cx)**2)
            phi[i, j] = np.arctan2((i - cy), (j - cx))   
#################################################################



#################################################################

def circle(r, a, Ny=10, Nx=10):
    inside = np.zeros((Ny, Nx), dtype=bool)
    for i in range(Ny):
        for j in range(Nx):
            if r[i, j] <= a:
                inside[i, j] = True
    return inside
#################################################################





#################################################################

def square(a, Ny=10, Nx=10, cy=0, cx=0):
    inside = np.zeros((Ny, Nx), dtype=bool)
    for i in range(Ny):
        for j in range(Nx):
            Y = i - cy
            X = j - cx
            if (X <= a and X >= -a and Y <= a and Y >= -a):
                inside[i, j] = True
    return inside
#################################################################



#################################################################

def hexagon(a, Ny=10, Nx=10, cy=0, cx=0):
    inside = np.zeros((Ny, Nx), dtype=bool)
    for i in range(Ny):
        for j in range(Nx):
            Y = i - cy
            X = j - cx
            if (Y <= X + 2*a and Y >= - X - 2*a and Y >= X - 2*a and Y <= -X + 2*a and X >= - a and X <= a):
                inside[i, j] = True
    return inside
#################################################################   




#################################################################

def corrugatedEllipse(r, phi, a, Ny=10, Nx=10, A=1, N=0, epsilon=0):
    inside = np.zeros((Ny, Nx), dtype=bool)
    for i in range(Ny):
        for j in range(Nx):
            if (r[i, j] <= a * (1 / ((np.cos(phi[i, j]))**2 + (A * np.sin(phi[i, j]))**2) + epsilon * np.cos(N * phi[i, j]))):
                inside[i, j] = True
    return inside
#################################################################

                


#################################################################

def JanusHalf(r, er, a, theta, er2, er3, Ny=10, Nx=10, cy=0, cx=0):
    for i in range(Ny):
        for j in range(Nx):
            Y = i - cy
            X = j - cx
            if (theta == 0 or theta == 360):
                if (r[i, j] <= a and X <= 0):
                    er[i, j]  = er2
                elif (r[i, j] <= a and X > 0):
                    er[i, j]  = er3
            elif (theta == 180):
                if (r[i, j] <= a and X >= 0):
                    er[i, j]  = er2
                elif (r[i, j] <= a and X < 0):
                    er[i, j]  = er3
            elif (0 < theta < 180):
                if (r[i, j] <= a and Y >= X * np.tan(np.pi/2 - theta * np.pi/180)):
                    er[i, j]  = er2
                elif (r[i, j] <= a and Y < X * np.tan(np.pi/2 - theta * np.pi/180)):
                    er[i, j]  = er3
            elif (180 < theta < 360):
                if (r[i, j] <= a and Y <= X * np.tan(3*np.pi/2 - theta * np.pi/180)):
                    er[i, j]  = er2
                elif (r[i, j] <= a and Y > X * np.tan(3*np.pi/2 - theta * np.pi/180)):
                    er[i, j]  = er3
#################################################################




#################################################################

def JanusSegment(r, er, a, er2, er3, Ny=10, Nx=10, cy=0, cx=0, shift=0):
    for i in range(Ny):
        for j in range(Nx):
            Y = i - cy
            X = j - cx
            if (r[i, j] <= a and X <= -(a-shift)):
                er[i, j]  = er2
            elif (r[i, j] <= a and X > -(a-shift)):
                er[i, j]  = er3
#################################################################




#################################################################

def JanusSector(r, er, a, er2, er3, Ny=10, Nx=10, cy=0, cx=0, theta=0):
    for i in range(Ny):
        for j in range(Nx):
            Y = i - cy
            X = j - cx
            if (r[i, j] <= a):
                if (theta <= 90):
                    if (X <= 0 and np.tan(theta * np.pi/180)*X <= Y <= -np.tan(theta * np.pi/180)*X):
                        er[i, j]  = er2
                    else :
                        er[i, j]  = er3
                else :
                    if (X > 0 and -np.tan((180-theta) * np.pi/180)*X <= Y <= np.tan((180-theta) * np.pi/180)*X):
                        er[i, j]  = er3
                    else :
                        er[i, j]  = er2
#################################################################














