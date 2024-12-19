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

ExI, EyI, EzI, HxI, HyI, HzI = [initialize_field(Nz, Ny, Nx) for _ in range(6)]
Ex, Ey, Ez, Hx, Hy, Hz = [initialize_field(Nz, Ny, Nx) for _ in range(6)]
Ex_scat, Ey_scat, Ez_scat, Hx_scat, Hy_scat, Hz_scat = [initialize_field(Nz, Ny, Nx) for _ in range(6)]


'''initializing the distribution functions of electric and magnetic fields'''
def initilize_dis_func(Nz=10, Ny=10, Nx=10, Q=7):
    return np.zeros((Nz, Ny, Nx, Q), dtype=np.float32, order='C')

exI, eyI, ezI, hxI, hyI, hzI = [initilize_dis_func(Nz, Ny, Nx, Q) for _ in range(6)]
exbI, eybI, ezbI, hxbI, hybI, hzbI = [initilize_dis_func(Nz, Ny, Nx, Q) for _ in range(6)]
ex, ey, ez, hx, hy, hz = [initilize_dis_func(Nz, Ny, Nx, Q) for _ in range(6)]
exb, eyb, ezb, hxb, hyb, hzb = [initilize_dis_func(Nz, Ny, Nx, Q) for _ in range(6)]

'''initializing the material properties for incident fields'''
erI, murI = initialize_material_properties(er1, mur1, Nz, Ny, Nx)
###############################################################################################################







t1 = time.time()


fig = plt.figure(figsize = (3.35*3, 3.35*2), dpi=600)
gs=GridSpec(2,3)
plt.ion()

fps = 10

Video_Name = "EzTotScat.mp4"

FFMpegWriter = animation.writers['ffmpeg']

metadata = dict(title='Ez', artist='', comment='Movie support!')
writer = FFMpegWriter(fps=fps, metadata=metadata)

with writer.saving(fig, Video_Name,300):
    for t in range(int(Time)):
        
        

        #################################################################################################################
        ########                                         LBM CALCULATION                                          #######
        #################################################################################################################

        '''initialization of macroscopic fields'''
        myclib.initializeField(ExI, EyI, EzI, HxI, HyI, HzI, Nz, Ny, Nx, N)
        myclib.initializeField(Ex, Ey, Ez, Hx, Hy, Hz, Nz, Ny, Nx, N)

        '''computation of macroscopic fields from distribution function'''
        myclib.macroField(exI, erI, ExI, Nz, Ny, Nx, Q, N)
        myclib.macroField(eyI, erI, EyI, Nz, Ny, Nx, Q, N)
        myclib.macroField(ezI, erI, EzI, Nz, Ny, Nx, Q, N)

        myclib.macroField(hxI, murI, HxI, Nz, Ny, Nx, Q, N)
        myclib.macroField(hyI, murI, HyI, Nz, Ny, Nx, Q, N)
        myclib.macroField(hzI, murI, HzI, Nz, Ny, Nx, Q, N)


        myclib.macroField(ex, er, Ex, Nz, Ny, Nx, Q, N)
        myclib.macroField(ey, er, Ey, Nz, Ny, Nx, Q, N)
        myclib.macroField(ez, er, Ez, Nz, Ny, Nx, Q, N)

        myclib.macroField(hx, mur, Hx, Nz, Ny, Nx, Q, N)
        myclib.macroField(hy, mur, Hy, Nz, Ny, Nx, Q, N)
        myclib.macroField(hz, mur, Hz, Nz, Ny, Nx, Q, N)


        if (t >= 0):
            
            '''source wave'''
            planeWaveTM(EzI, HyI, t, omega, xloc, ymin, ymax, zmin, zmax)
            planeWaveTM(Ez, Hy, t, omega, xloc, ymin, ymax, zmin, zmax)

            '''calculation of scattered fields'''
            Ez_scat = Ez - EzI
            Hx_scat = Hx - HxI
            Hy_scat = Hy - HyI

##            Ez_scat[scatterer] = 0
##            Hx_scat[scatterer] = 0
##            Hy_scat[scatterer] = 0

            '''collision and streaming (the 2 steps of LBM) when field is forced'''
            myclib.collForcingNode(exI, eyI, ezI, hxI, hyI, hzI, exbI, eybI, ezbI, hxbI, hybI, hzbI, ExI, EyI, EzI, HxI, HyI, HzI, erI, murI, Nz, Ny, Nx, Q, xloc, ymin, ymax, zmin, zmax, N)
            myclib.collForcingNode(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ex, Ey, Ez, Hx, Hy, Hz, er, mur, Nz, Ny, Nx, Q, xloc, ymin, ymax, zmin, zmax, N)
            myclib.streaming(exI, eyI, ezI, hxI, hyI, hzI, exbI, eybI, ezbI, hxbI, hybI, hzbI, Nz, Ny, Nx, Q, N)
            myclib.streaming(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Nz, Ny, Nx, Q, N)

        else:     
            '''collision and streaming (the 2 steps of LBM) when field is not forced'''
            myclib.collision(exI, eyI, ezI, hxI, hyI, hzI, exbI, eybI, ezbI, hxbI, hybI, hzbI, ExI, EyI, EzI, HxI, HyI, HzI, erI, murI, Nz, Ny, Nx, Q, N)
            myclib.collision(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ex, Ey, Ez, Hx, Hy, Hz, er, mur, Nz, Ny, Nx, Q, N)
            myclib.streaming(exI, eyI, ezI, hxI, hyI, hzI, exbI, eybI, ezbI, hxbI, hybI, hzbI, Nz, Ny, Nx, Q, N)
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





            ax4 = fig.add_subplot(gs[1,0])
            plt.title(r'$E_z^{scat}$')
            im4 = plt.imshow(Ez_scat[Nz//2, :, :], vmin = -1, vmax = 1, cmap='seismic', origin='lower')
            im = plt.imshow(er[Nz//2, :, :], extent=(0, Nx, 0, Ny), cmap='binary', origin='lower', alpha=0.1)
            
            ax4.set_xticks(np.linspace(0,Nx,4))
            ax4.set_xticklabels([])
            ax4.set_xlabel(r'$x$', fontsize=12)
            ax4.set_yticks(np.linspace(0,Ny,4))
            ax4.set_yticklabels([])
            ax4.set_ylabel(r'$y$', fontsize=12)


            ax5 = fig.add_subplot(gs[1,1])
            plt.title(r'$E_z^{scat}$')
            im5 = plt.imshow(Ez_scat[:, :, Nx//2], vmin = -1, vmax = 1, cmap='seismic', origin='lower')
            im = plt.imshow(er[:, :, Nx//2], extent=(0, Nx, 0, Ny), cmap='binary', origin='lower', alpha=0.1)
            
            ax5.set_xticks(np.linspace(0,Nx,4))
            ax5.set_xticklabels([])
            ax5.set_xlabel(r'$y$', fontsize=12)
            ax5.set_yticks(np.linspace(0,Ny,4))
            ax5.set_yticklabels([])
            ax5.set_ylabel(r'$z$', fontsize=12)



            ax6 = fig.add_subplot(gs[1,2])
            plt.title(r'$E_z^{scat}$')
            im6 = plt.imshow(Ez_scat[:, Ny//2, :], vmin = -1, vmax = 1, cmap='seismic', origin='lower')
            im = plt.imshow(er[:, Ny//2, :], extent=(0, Nx, 0, Ny), cmap='binary', origin='lower', alpha=0.1)
            
            ax6.set_xticks(np.linspace(0,Nx,4))
            ax6.set_xticklabels([])
            ax6.set_xlabel(r'$x$', fontsize=12)
            ax6.set_yticks(np.linspace(0,Ny,4))
            ax6.set_yticklabels([])
            ax6.set_ylabel(r'$z$', fontsize=12)








##            plt.colorbar(im2, location='right', shrink=0.97, aspect=20)
       
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
