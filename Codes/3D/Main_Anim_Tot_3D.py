import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec

from Module_Traction_3D import *
from Module_EM_Wave_3D import *
from Module_Parameters_3D import *
from Module_Shared_Lib_3D import *


t0 = time.time()


pictures = 'gallery'
if not os.path.exists(pictures):
    os.makedirs(pictures)




print("Number of time steps :", int(Time))




###############################################################################################################
####      DEFINING AND INITILIZING VARIABLES FOR MACROSCOPIC FIELDS AND DISTRIBUTION FUNCTIONS             ####
###############################################################################################################

'''initializing the electric and magnetic fields'''
def initialize_field(Nz=10, Ny=10, Nx=10):
    return np.zeros((Nz, Ny, Nx), dtype=np.float32, order='C')

Ex, Ey, Ez, Hx, Hy, Hz = [initialize_field(Nz, Ny, Nx) for _ in range(6)]


'''initializing the distribution functions of electric and magnetic fields'''
def initilize_dis_func(Nz=10, Ny=10, Nx=10, Q=7):
    return np.zeros((Nz, Ny, Nx, Q), dtype=np.float32, order='C')

ex, ey, ez, hx, hy, hz = [initilize_dis_func(Nz, Ny, Nx, Q) for _ in range(6)]
exb, eyb, ezb, hxb, hyb, hzb = [initilize_dis_func(Nz, Ny, Nx, Q) for _ in range(6)]

###############################################################################################################









t1 = time.time()


fig = plt.figure(figsize = (3.35*3, 3.35*1), dpi=600)
gs=GridSpec(1,3)
plt.ion()

fps = 10

Video_Name = "Ez_tot.mp4"

FFMpegWriter = animation.writers['ffmpeg']

metadata = dict(title='Ez_tot', artist='', comment='Movie support!')
writer = FFMpegWriter(fps=fps, metadata=metadata)

with writer.saving(fig, Video_Name,300):
    for t in range(int(Time)):
        
        

        #################################################################################################################
        ########                                         LBM CALCULATION                                          #######
        #################################################################################################################

        '''initialization of macroscopic fields'''
        myclib.initializeField(Ex, Ey, Ez, Hx, Hy, Hz, Nz, Ny, Nx, N)

        '''computation of macroscopic fields from distribution function'''
        myclib.macroField(ex, er, Ex, Nz, Ny, Nx, Q, N)
        myclib.macroField(ey, er, Ey, Nz, Ny, Nx, Q, N)
        myclib.macroField(ez, er, Ez, Nz, Ny, Nx, Q, N)

        myclib.macroField(hx, mur, Hx, Nz, Ny, Nx, Q, N)
        myclib.macroField(hy, mur, Hy, Nz, Ny, Nx, Q, N)
        myclib.macroField(hz, mur, Hz, Nz, Ny, Nx, Q, N)

        
        if (t >= 0):
            
            '''source wave'''
            planeWaveTM(Ez, Hy, t, omega, xloc, ymin, ymax, zmin, zmax)

            '''collision and streaming (the 2 steps of LBM) when field is forced'''
            myclib.collForcingNode(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ex, Ey, Ez, Hx, Hy, Hz, er, mur, Nz, Ny, Nx, Q, xloc, ymin, ymax, zmin, zmax, N)
            myclib.streaming(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Nz, Ny, Nx, Q, N)
        else:    
            '''collision and streaming (the 2 steps of LBM) when field is not forced'''
            myclib.collision(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ex, Ey, Ez, Hx, Hy, Hz, er, mur, Nz, Ny, Nx, Q, N)
            myclib.streaming(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Nz, Ny, Nx, Q, N)
                
        ###############################################################################################################


        

    ###############################################################################################################
    ##########                                           ANIMATION                                       ##########
    ###############################################################################################################

        if (t%10 == 0):
                
            fig.clear()

            ########### Plot parameters ################
            plt.rc('font', family = 'serif', size = 10)
            plt.rc('xtick', labelsize = 10)
            plt.rc('ytick', labelsize = 10)
            plt.rc('lines', markersize = 2, lw = 0.75)
            plt.rc('text', usetex = True)
            plt.rcParams['image.cmap']='seismic'
            ###########################################
            
            
                        
            ax1 = fig.add_subplot(gs[0,0])
            plt.title(r'$E_z^{tot}$')
            im1 = plt.imshow(Ez[Nz//2, :, :], vmin = -1, vmax = 1, cmap='seismic', origin='lower')
            im = plt.imshow(er[Nz//2, :, :], extent=(0, Nx, 0, Ny), cmap='binary', origin='lower', alpha=0.1)
            
            ax1.set_xticks(np.linspace(0,Nx,4))
            ax1.set_xticklabels([])
            ax1.set_xlabel(r'$x$', fontsize=12)
            ax1.set_yticks(np.linspace(0,Ny,4))
            ax1.set_yticklabels([])
            ax1.set_ylabel(r'$y$', fontsize=12)


            ax2 = fig.add_subplot(gs[0,1])
            plt.title(r'$E_z^{tot}$')
            im2 = plt.imshow(Ez[:, :, Nx//2], vmin = -1, vmax = 1, cmap='seismic', origin='lower')
            im = plt.imshow(er[:, :, Nx//2], extent=(0, Nx, 0, Ny), cmap='binary', origin='lower', alpha=0.1)
            
            ax2.set_xticks(np.linspace(0,Nx,4))
            ax2.set_xticklabels([])
            ax2.set_xlabel(r'$y$', fontsize=12)
            ax2.set_yticks(np.linspace(0,Ny,4))
            ax2.set_yticklabels([])
            ax2.set_ylabel(r'$z$', fontsize=12)



            ax3 = fig.add_subplot(gs[0,2])
            plt.title(r'$E_z^{tot}$')
            im3 = plt.imshow(Ez[:, Ny//2, :], vmin = -1, vmax = 1, cmap='seismic', origin='lower')
            im = plt.imshow(er[:, Ny//2, :], extent=(0, Nx, 0, Ny), cmap='binary', origin='lower', alpha=0.1)
            
            ax3.set_xticks(np.linspace(0,Nx,4))
            ax3.set_xticklabels([])
            ax3.set_xlabel(r'$x$', fontsize=12)
            ax3.set_yticks(np.linspace(0,Ny,4))
            ax3.set_yticklabels([])
            ax3.set_ylabel(r'$z$', fontsize=12)


##            plt.colorbar(im1, location='right', shrink=0.97, aspect=20)
       
##            plt.savefig(pictures+"/pic."+str(t).zfill(4)+".png")
            writer.grab_frame()




    ###############################################################################################################
    ###############################################################################################################



###############################################################################################################
        t2 = time.time()
        
        if (t > 0 and t%100 == 0):
            remaining_time = (t2 - t1) * (int(Time) - t) / (t*60)
            print(f"Approximate time left: {remaining_time:.2f} minutes", end="\r")

t3 = time.time()
total_time = (t3 - t0) / 60
print(f"\nTotal time taken: {total_time:.2f} minutes")
###############################################################################################################
