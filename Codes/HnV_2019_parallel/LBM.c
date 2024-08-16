#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<omp.h>


#define W 1.0/4.0

// number of parallel threads
int N;




/*---------------------- SET OF VECTORS ------------------------*/

// velocity vector
int vx[12] = { 1, -1,  0,  0,  0,  0,     -1,  1,  0,  0,  0,  0};
int vy[12] = { 0,  0,  1, -1,  0,  0,      0,  0, -1,  1,  0,  0};
int vz[12] = { 0,  0,  0,  0,  1, -1,      0,  0,  0,  0, -1,  1};

// construction vector of electric field
int ex[12] = { 0,  0,  0,  1,  0, -1,      0,  0,  0,  1,  0, -1};
int ey[12] = { 0,  1,  0,  0, -1,  0,      0,  1,  0,  0, -1,  0};
int ez[12] = { 1,  0, -1,  0,  0,  0,      1,  0, -1,  0,  0,  0};

// construction vector of magnetic field
int hx[12] = { 0,  0, -1,  0,  1,  0,      0,  0,  1,  0, -1,  0};
int hy[12] = {-1,  0,  0,  0,  0,  1,      1,  0,  0,  0,  0, -1};
int hz[12] = { 0, -1,  0,  1,  0,  0,      0,  1,  0, -1,  0,  0};

/*--------------------------------------------------------------*/


/*---------------------------------------------------------------------------------------------*/
// calculates the macroscopic fields from the distribution functions


void macroField(float* f, float* Ex, float* Ey, float* Ez, float* Hx, float* Hy, float* Hz, 
					float* Px, float* Py, float* Pz, float* Mx, float* My, float* Mz, float* er, float* mur, int y, int x, int q, int N) {
	
	int i, j, k;
	
	for(i = 0; i < y; i++) {
	
		omp_set_num_threads(N);
		#pragma omp parallel for private(k)
		
		for(j = 0; j < x; j++) {
		
			float sum_Ex = 0, sum_Ey = 0, sum_Ez = 0;
			float sum_Hx = 0, sum_Hy = 0, sum_Hz = 0;
			
			for(k = 0; k < q; k++) {
				
				// summing of the distribution functions over the velocity set
				
				sum_Ex += f[i*x*q + j*q + k] * ex[k];
				sum_Ey += f[i*x*q + j*q + k] * ey[k];
				sum_Ez += f[i*x*q + j*q + k] * ez[k];
				
				sum_Hx += f[i*x*q + j*q + k] * hx[k];
				sum_Hy += f[i*x*q + j*q + k] * hy[k];
				sum_Hz += f[i*x*q + j*q + k] * hz[k];
			}
			
			// macroscopic field calculation
			
			Ex[i*x + j] = (sum_Ex + Px[i*x + j]) / er[i*x + j];
			Ey[i*x + j] = (sum_Ey + Py[i*x + j]) / er[i*x + j];
			Ez[i*x + j] = (sum_Ez + Pz[i*x + j]) / er[i*x + j];
			
			Hx[i*x + j] = (sum_Hx + Mx[i*x + j]) / mur[i*x + j];
			Hy[i*x + j] = (sum_Hy + My[i*x + j]) / mur[i*x + j];
			Hz[i*x + j] = (sum_Hz + Mz[i*x + j]) / mur[i*x + j];
		}
	}
}
/*--------------------------------------------------------------------------------------------*/






/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------COLLISION IF FIELDS ARE FORCED------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/


void collForcingNode(float* f, float* fb, float* Ex, float* Ey, float* Ez, float* Hx, float* Hy, float* Hz, 
					float* Px, float* Py, float* Pz, float* Mx, float* My, float* Mz, float* Pxb, float* Pyb, float* Pzb, float* Mxb, float* Myb, float* Mzb, 
					float* er, float* mur, int y, int x, int q, int xloc, int ymin, int ymax, int N) {

	int i, j;
	
	float feq0, feq1, feq2, feq3, feq4, feq5, feq6, feq7, feq8, feq9, feq10, feq11;
	float PxEq, PyEq, PzEq, MxEq, MyEq, MzEq;

	for (i = 0; i < y; i++) {
	
		omp_set_num_threads(N);
		#pragma omp parallel for
		
		for (j = 0; j < x; j++) {
			
			int index   = i*x + j;
			int index_f = i*x*q + j*q;

			// equilibrium distribution functions
			feq0  = W * (Ex[index]*ex[0]  + Hx[index]*hx[0]  + Ey[index]*ey[0]  + Hy[index]*hy[0]  + Ez[index]*ez[0]  + Hz[index]*hz[0]);
			feq1  = W * (Ex[index]*ex[1]  + Hx[index]*hx[1]  + Ey[index]*ey[1]  + Hy[index]*hy[1]  + Ez[index]*ez[1]  + Hz[index]*hz[1]);
			feq2  = W * (Ex[index]*ex[2]  + Hx[index]*hx[2]  + Ey[index]*ey[2]  + Hy[index]*hy[2]  + Ez[index]*ez[2]  + Hz[index]*hz[2]);
			feq3  = W * (Ex[index]*ex[3]  + Hx[index]*hx[3]  + Ey[index]*ey[3]  + Hy[index]*hy[3]  + Ez[index]*ez[3]  + Hz[index]*hz[3]);
			feq4  = W * (Ex[index]*ex[4]  + Hx[index]*hx[4]  + Ey[index]*ey[4]  + Hy[index]*hy[4]  + Ez[index]*ez[4]  + Hz[index]*hz[4]);
			feq5  = W * (Ex[index]*ex[5]  + Hx[index]*hx[5]  + Ey[index]*ey[5]  + Hy[index]*hy[5]  + Ez[index]*ez[5]  + Hz[index]*hz[5]);
			feq6  = W * (Ex[index]*ex[6]  + Hx[index]*hx[6]  + Ey[index]*ey[6]  + Hy[index]*hy[6]  + Ez[index]*ez[6]  + Hz[index]*hz[6]);
			feq7  = W * (Ex[index]*ex[7]  + Hx[index]*hx[7]  + Ey[index]*ey[7]  + Hy[index]*hy[7]  + Ez[index]*ez[7]  + Hz[index]*hz[7]);
			feq8  = W * (Ex[index]*ex[8]  + Hx[index]*hx[8]  + Ey[index]*ey[8]  + Hy[index]*hy[8]  + Ez[index]*ez[8]  + Hz[index]*hz[8]);
			feq9  = W * (Ex[index]*ex[9]  + Hx[index]*hx[9]  + Ey[index]*ey[9]  + Hy[index]*hy[9]  + Ez[index]*ez[9]  + Hz[index]*hz[9]);
			feq10 = W * (Ex[index]*ex[10] + Hx[index]*hx[10] + Ey[index]*ey[10] + Hy[index]*hy[10] + Ez[index]*ez[10] + Hz[index]*hz[10]);
			feq11 = W * (Ex[index]*ex[11] + Hx[index]*hx[11] + Ey[index]*ey[11] + Hy[index]*hy[11] + Ez[index]*ez[11] + Hz[index]*hz[11]);
			
			// equilibrium polarization
			PxEq = (er[index] - 1) * Ex[index];
			PyEq = (er[index] - 1) * Ey[index];
			PzEq = (er[index] - 1) * Ez[index];
			
			// equilibrium magnetization
			MxEq = (mur[index] - 1) * Hx[index];
			MyEq = (mur[index] - 1) * Hy[index];
			MzEq = (mur[index] - 1) * Hz[index];

			
			
			if (i >= ymin && i < ymax && j == xloc) { // no collision at the forcing node
				
				// polarization at forcing nodes
				Pxb[index] = PxEq;
				Pyb[index] = PyEq;
				Pzb[index] = PzEq;
				
				// magnetization at forcing nodes
				Mxb[index] = MxEq;
				Myb[index] = MyEq;
				Mzb[index] = MzEq;
				
				// distribution functions at forcing nodes
				fb[index_f + 0]  = feq0;	
				fb[index_f + 1]  = feq1;
				fb[index_f + 2]  = feq2;
				fb[index_f + 3]  = feq3;
				fb[index_f + 4]  = feq4;
				fb[index_f + 5]  = feq5;
				fb[index_f + 6]  = feq6;	
				fb[index_f + 7]  = feq7;
				fb[index_f + 8]  = feq8;
				fb[index_f + 9]  = feq9;
				fb[index_f + 10] = feq10;
				fb[index_f + 11] = feq11;
			}
			

			else {
				
				// polarization at non-forcing nodes
				Pxb[index] = 2 * PxEq - Px[index];
				Pyb[index] = 2 * PyEq - Py[index];
				Pzb[index] = 2 * PzEq - Pz[index];
				
				// magnetization at non-forcing nodes
				Mxb[index] = 2 * MxEq - Mx[index];
				Myb[index] = 2 * MyEq - My[index];
				Mzb[index] = 2 * MzEq - Mz[index];

				// distribution functions at non-forcing nodes
				fb[index_f + 0]  = 2 * feq0  - f[index_f + 0];
				fb[index_f + 1]  = 2 * feq1  - f[index_f + 1];
				fb[index_f + 2]  = 2 * feq2  - f[index_f + 2];
				fb[index_f + 3]  = 2 * feq3  - f[index_f + 3];
				fb[index_f + 4]  = 2 * feq4  - f[index_f + 4];
				fb[index_f + 5]  = 2 * feq5  - f[index_f + 5];
				fb[index_f + 6]  = 2 * feq6  - f[index_f + 6];
				fb[index_f + 7]  = 2 * feq7  - f[index_f + 7];
				fb[index_f + 8]  = 2 * feq8  - f[index_f + 8];
				fb[index_f + 9]  = 2 * feq9  - f[index_f + 9];
				fb[index_f + 10] = 2 * feq10 - f[index_f + 10];
				fb[index_f + 11] = 2 * feq11 - f[index_f + 11];
			}
		}
	}
}




/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/










/*---------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------COLLISION IF FIELDS ARE NOT FORCED----------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/


void collNotForcingNode(float* f, float* fb, float* Ex, float* Ey, float* Ez, float* Hx, float* Hy, float* Hz, float* Px, float* Py, float* Pz, 
						float* Mx, float* My, float* Mz, float* Pxb, float* Pyb, float* Pzb, float* Mxb, float* Myb, float* Mzb, float* er, float* mur, 
						int y, int x, int q, int N) {

	int i, j;
	
	float feq0, feq1, feq2, feq3, feq4, feq5, feq6, feq7, feq8, feq9, feq10, feq11;
	float PxEq, PyEq, PzEq, MxEq, MyEq, MzEq;	

	for (i = 0; i < y; i++) {
		omp_set_num_threads(N);
		#pragma omp parallel for
		for (j = 0; j < x; j++) {
			

			int index   = i*x + j;
			int index_f = i*x*q + j*q;

			// equilibrium distribution functions
			feq0  = W * (Ex[index]*ex[0]  + Hx[index]*hx[0]  + Ey[index]*ey[0]  + Hy[index]*hy[0]  + Ez[index]*ez[0]  + Hz[index]*hz[0]);
			feq1  = W * (Ex[index]*ex[1]  + Hx[index]*hx[1]  + Ey[index]*ey[1]  + Hy[index]*hy[1]  + Ez[index]*ez[1]  + Hz[index]*hz[1]);
			feq2  = W * (Ex[index]*ex[2]  + Hx[index]*hx[2]  + Ey[index]*ey[2]  + Hy[index]*hy[2]  + Ez[index]*ez[2]  + Hz[index]*hz[2]);
			feq3  = W * (Ex[index]*ex[3]  + Hx[index]*hx[3]  + Ey[index]*ey[3]  + Hy[index]*hy[3]  + Ez[index]*ez[3]  + Hz[index]*hz[3]);
			feq4  = W * (Ex[index]*ex[4]  + Hx[index]*hx[4]  + Ey[index]*ey[4]  + Hy[index]*hy[4]  + Ez[index]*ez[4]  + Hz[index]*hz[4]);
			feq5  = W * (Ex[index]*ex[5]  + Hx[index]*hx[5]  + Ey[index]*ey[5]  + Hy[index]*hy[5]  + Ez[index]*ez[5]  + Hz[index]*hz[5]);
			feq6  = W * (Ex[index]*ex[6]  + Hx[index]*hx[6]  + Ey[index]*ey[6]  + Hy[index]*hy[6]  + Ez[index]*ez[6]  + Hz[index]*hz[6]);
			feq7  = W * (Ex[index]*ex[7]  + Hx[index]*hx[7]  + Ey[index]*ey[7]  + Hy[index]*hy[7]  + Ez[index]*ez[7]  + Hz[index]*hz[7]);
			feq8  = W * (Ex[index]*ex[8]  + Hx[index]*hx[8]  + Ey[index]*ey[8]  + Hy[index]*hy[8]  + Ez[index]*ez[8]  + Hz[index]*hz[8]);
			feq9  = W * (Ex[index]*ex[9]  + Hx[index]*hx[9]  + Ey[index]*ey[9]  + Hy[index]*hy[9]  + Ez[index]*ez[9]  + Hz[index]*hz[9]);
			feq10 = W * (Ex[index]*ex[10] + Hx[index]*hx[10] + Ey[index]*ey[10] + Hy[index]*hy[10] + Ez[index]*ez[10] + Hz[index]*hz[10]);
			feq11 = W * (Ex[index]*ex[11] + Hx[index]*hx[11] + Ey[index]*ey[11] + Hy[index]*hy[11] + Ez[index]*ez[11] + Hz[index]*hz[11]);
			
			// equilibrium polarization
			PxEq = (er[index] - 1) * Ex[index];
			PyEq = (er[index] - 1) * Ey[index];
			PzEq = (er[index] - 1) * Ez[index];
			
			// equilibrium magnetization
			MxEq = (mur[index] - 1) * Hx[index];
			MyEq = (mur[index] - 1) * Hy[index];
			MzEq = (mur[index] - 1) * Hz[index];

			// polarization at non-forcing nodes
			Pxb[index] = 2 * PxEq - Px[index];
			Pyb[index] = 2 * PyEq - Py[index];
			Pzb[index] = 2 * PzEq - Pz[index];
			
			// magnetization at non-forcing nodes
			Mxb[index] = 2 * MxEq - Mx[index];
			Myb[index] = 2 * MyEq - My[index];
			Mzb[index] = 2 * MzEq - Mz[index];
			
			// distribution functions at non-forcing nodes
			fb[index_f + 0]  = 2 * feq0  - f[index_f + 0];
			fb[index_f + 1]  = 2 * feq1  - f[index_f + 1];
			fb[index_f + 2]  = 2 * feq2  - f[index_f + 2];
			fb[index_f + 3]  = 2 * feq3  - f[index_f + 3];
			fb[index_f + 4]  = 2 * feq4  - f[index_f + 4];
			fb[index_f + 5]  = 2 * feq5  - f[index_f + 5];
			fb[index_f + 6]  = 2 * feq6  - f[index_f + 6];
			fb[index_f + 7]  = 2 * feq7  - f[index_f + 7];
			fb[index_f + 8]  = 2 * feq8  - f[index_f + 8];
			fb[index_f + 9]  = 2 * feq9  - f[index_f + 9];
			fb[index_f + 10] = 2 * feq10 - f[index_f + 10];
			fb[index_f + 11] = 2 * feq11 - f[index_f + 11];
		}
	}
}


/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/







/*---------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------- STREAMING ---------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/


void streaming(float* f, float* fb, float* Px, float* Py, float* Pz, float* Mx, float* My, float* Mz, 
					float* Pxb, float* Pyb, float* Pzb, float* Mxb, float* Myb, float* Mzb, int y, int x, int q, int N) {

	int i, j;
	
	

	for (i = 0; i < y; i++) {
	
		omp_set_num_threads(N);
		#pragma omp parallel for
		
		for (j = 0; j < x; j++) {
				
			int index   = i*x + j;
			int index_f = i*x*q + j*q;
			
			
			int j1 = j - 1, j2 = j + 1;
			int i3 = i - 1, i4 = i + 1;
			
			
			if (j == 0)
				j1 = j;
			if (j == x-1)
				j2 = j;
			if (i == 0)
				i3 = i;
			if (i == y-1)
				i4 = i;

			
			
			Px[index] = Pxb[index];
			Py[index] = Pyb[index];
			Pz[index] = Pzb[index];
			
			Mx[index] = Mxb[index];
			My[index] = Myb[index];
			Mz[index] = Mzb[index];
			
			
			f[index_f + 0]  = fb[i*x*q + j1*q + 0];
			f[index_f + 1]  = fb[i*x*q + j2*q + 1];
			f[index_f + 2]  = fb[i3*x*q + j*q + 2];
			f[index_f + 3]  = fb[i4*x*q + j*q + 3];
			f[index_f + 4]  = fb[index_f  + 4];
			f[index_f + 5]  = fb[index_f  + 5];
			f[index_f + 6]  = fb[i*x*q + j2*q + 6];
			f[index_f + 7]  = fb[i*x*q + j1*q + 7];
			f[index_f + 8]  = fb[i4*x*q + j*q + 8];
			f[index_f + 9]  = fb[i3*x*q + j*q + 9];
			f[index_f + 10] = fb[index_f + 10];
			f[index_f + 11] = fb[index_f + 11];
		}
	}
}


/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/




