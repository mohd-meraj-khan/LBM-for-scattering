#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<omp.h>


#define W 1.0/6.0

// number of threads in parallel
int N = 2;


int V[21] = {0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1};


/*---------------------------------------------------------------------------------------------*/
// calculates the macroscopic fields from the distribution functions


void macroField(float* dis_func, float* mat_prop, float* field, int y, int x, int q) {
	
	int i, j, k;
	
	
	for(i = 0; i < y; i++) {
		omp_set_num_threads(N);
		#pragma omp parallel for private(k)
		for(j = 0; j < x; j++) {
			for(k = 0; k < q; k++) {
				field[i*x + j] += dis_func[i*x*q + j*q + k] / mat_prop[i*x + j];
			}
		}
	}
}
/*--------------------------------------------------------------------------------------------*/






/*----------------------------------------------------------------------------------------------------*/
// inintialize the macroscopic field


void initializeField(float* Ex, float* Ey, float* Ez, float* Hx, float* Hy, float* Hz, int y, int x) {
	
	int i, j;

	for(i = 0; i < y; i++){
		for(j = 0; j < x; j++){				
			Ex[i*x + j] = 0;
			Ey[i*x + j] = 0;
			Ez[i*x + j] = 0;
			
			Hx[i*x + j] = 0;
			Hy[i*x + j] = 0;
			Hz[i*x + j] = 0;
		}
	}
}

/*----------------------------------------------------------------------------------------------------*/








/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------COLLISION + STREAMING IF FIELDS ARE FORCED------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/


void collForcingNode(float* ex, float* ey, float* ez, float* hx, float* hy, float* hz, float* exb, float* eyb, float* ezb, float* hxb, float* hyb, float* hzb, float* Ex, float* Ey, float* Ez, float* Hx, float* Hy, float* Hz, float* er, float* mur, int y, int x, int q, int xloc, int ymin, int ymax) {

	int i, j;
	
	
	for (i = 0; i < y; i++) {
		omp_set_num_threads(N);
		#pragma omp parallel for
		for (j = 0; j < x; j++) {

			
			if (i >= ymin && i < ymax && j == xloc) { // no collision at the forcing node

				exb[i*x*q + j*q + 0] = ((er[i*x + j] - 1) * Ex[i*x + j]);	
				exb[i*x*q + j*q + 1] = W * (Ex[i*x + j] - V[1*3 + 1] * Hz[i*x + j] + V[1*3 + 2] * Hy[i*x + j]);
				exb[i*x*q + j*q + 2] = W * (Ex[i*x + j] - V[2*3 + 1] * Hz[i*x + j] + V[2*3 + 2] * Hy[i*x + j]);
				exb[i*x*q + j*q + 3] = W * (Ex[i*x + j] - V[3*3 + 1] * Hz[i*x + j] + V[3*3 + 2] * Hy[i*x + j]);
				exb[i*x*q + j*q + 4] = W * (Ex[i*x + j] - V[4*3 + 1] * Hz[i*x + j] + V[4*3 + 2] * Hy[i*x + j]);
				exb[i*x*q + j*q + 5] = W * (Ex[i*x + j] - V[5*3 + 1] * Hz[i*x + j] + V[5*3 + 2] * Hy[i*x + j]);
				exb[i*x*q + j*q + 6] = W * (Ex[i*x + j] - V[6*3 + 1] * Hz[i*x + j] + V[6*3 + 2] * Hy[i*x + j]);
				
				
				eyb[i*x*q + j*q + 0] = ((er[i*x + j] - 1) * Ey[i*x + j]);
				eyb[i*x*q + j*q + 1] = W * (Ey[i*x + j] - V[1*3 + 2] * Hx[i*x + j] + V[1*3 + 0] * Hz[i*x + j]);
				eyb[i*x*q + j*q + 2] = W * (Ey[i*x + j] - V[2*3 + 2] * Hx[i*x + j] + V[2*3 + 0] * Hz[i*x + j]);
				eyb[i*x*q + j*q + 3] = W * (Ey[i*x + j] - V[3*3 + 2] * Hx[i*x + j] + V[3*3 + 0] * Hz[i*x + j]);
				eyb[i*x*q + j*q + 4] = W * (Ey[i*x + j] - V[4*3 + 2] * Hx[i*x + j] + V[4*3 + 0] * Hz[i*x + j]);
				eyb[i*x*q + j*q + 5] = W * (Ey[i*x + j] - V[5*3 + 2] * Hx[i*x + j] + V[5*3 + 0] * Hz[i*x + j]);
				eyb[i*x*q + j*q + 6] = W * (Ey[i*x + j] - V[6*3 + 2] * Hx[i*x + j] + V[6*3 + 0] * Hz[i*x + j]);


				ezb[i*x*q + j*q + 0] = ((er[i*x + j] - 1) * Ez[i*x + j]);
				ezb[i*x*q + j*q + 1] = W * (Ez[i*x + j] - V[1*3 + 0] * Hy[i*x + j] + V[1*3 + 1] * Hx[i*x + j]);
				ezb[i*x*q + j*q + 2] = W * (Ez[i*x + j] - V[2*3 + 0] * Hy[i*x + j] + V[2*3 + 1] * Hx[i*x + j]);
				ezb[i*x*q + j*q + 3] = W * (Ez[i*x + j] - V[3*3 + 0] * Hy[i*x + j] + V[3*3 + 1] * Hx[i*x + j]);
				ezb[i*x*q + j*q + 4] = W * (Ez[i*x + j] - V[4*3 + 0] * Hy[i*x + j] + V[4*3 + 1] * Hx[i*x + j]);
				ezb[i*x*q + j*q + 5] = W * (Ez[i*x + j] - V[5*3 + 0] * Hy[i*x + j] + V[5*3 + 1] * Hx[i*x + j]);
				ezb[i*x*q + j*q + 6] = W * (Ez[i*x + j] - V[6*3 + 0] * Hy[i*x + j] + V[6*3 + 1] * Hx[i*x + j]);



				hxb[i*x*q + j*q + 0] = ((mur[i*x + j] - 1) * Hx[i*x + j]);
				hxb[i*x*q + j*q + 1] = W * (Hx[i*x + j] - V[1*3 + 2] * Ey[i*x + j] + V[1*3 + 1] * Ez[i*x + j]);
				hxb[i*x*q + j*q + 2] = W * (Hx[i*x + j] - V[2*3 + 2] * Ey[i*x + j] + V[2*3 + 1] * Ez[i*x + j]);
				hxb[i*x*q + j*q + 3] = W * (Hx[i*x + j] - V[3*3 + 2] * Ey[i*x + j] + V[3*3 + 1] * Ez[i*x + j]);
				hxb[i*x*q + j*q + 4] = W * (Hx[i*x + j] - V[4*3 + 2] * Ey[i*x + j] + V[4*3 + 1] * Ez[i*x + j]);
				hxb[i*x*q + j*q + 5] = W * (Hx[i*x + j] - V[5*3 + 2] * Ey[i*x + j] + V[5*3 + 1] * Ez[i*x + j]);
				hxb[i*x*q + j*q + 6] = W * (Hx[i*x + j] - V[6*3 + 2] * Ey[i*x + j] + V[6*3 + 1] * Ez[i*x + j]);
				
				
			
				hyb[i*x*q + j*q + 0] = ((mur[i*x + j] - 1) * Hy[i*x + j]);
				hyb[i*x*q + j*q + 1] = W * (Hy[i*x + j] - V[1*3 + 0] * Ez[i*x + j] + V[1*3 + 2] * Ex[i*x + j]);
				hyb[i*x*q + j*q + 2] = W * (Hy[i*x + j] - V[2*3 + 0] * Ez[i*x + j] + V[2*3 + 2] * Ex[i*x + j]);
				hyb[i*x*q + j*q + 3] = W * (Hy[i*x + j] - V[3*3 + 0] * Ez[i*x + j] + V[3*3 + 2] * Ex[i*x + j]);
				hyb[i*x*q + j*q + 4] = W * (Hy[i*x + j] - V[4*3 + 0] * Ez[i*x + j] + V[4*3 + 2] * Ex[i*x + j]);
				hyb[i*x*q + j*q + 5] = W * (Hy[i*x + j] - V[5*3 + 0] * Ez[i*x + j] + V[5*3 + 2] * Ex[i*x + j]);
				hyb[i*x*q + j*q + 6] = W * (Hy[i*x + j] - V[6*3 + 0] * Ez[i*x + j] + V[6*3 + 2] * Ex[i*x + j]);


				hzb[i*x*q + j*q + 0] = ((mur[i*x + j] - 1) * Hz[i*x + j]);
				hzb[i*x*q + j*q + 1] = W * (Hz[i*x + j] - V[1*3 + 1] * Ex[i*x + j] + V[1*3 + 0] * Ey[i*x + j]);
				hzb[i*x*q + j*q + 2] = W * (Hz[i*x + j] - V[2*3 + 1] * Ex[i*x + j] + V[2*3 + 0] * Ey[i*x + j]);
				hzb[i*x*q + j*q + 3] = W * (Hz[i*x + j] - V[3*3 + 1] * Ex[i*x + j] + V[3*3 + 0] * Ey[i*x + j]);
				hzb[i*x*q + j*q + 4] = W * (Hz[i*x + j] - V[4*3 + 1] * Ex[i*x + j] + V[4*3 + 0] * Ey[i*x + j]);
				hzb[i*x*q + j*q + 5] = W * (Hz[i*x + j] - V[5*3 + 1] * Ex[i*x + j] + V[5*3 + 0] * Ey[i*x + j]);
				hzb[i*x*q + j*q + 6] = W * (Hz[i*x + j] - V[6*3 + 1] * Ex[i*x + j] + V[6*3 + 0] * Ey[i*x + j]);
			
			}
			
							
			
			else {
				
				exb[i*x*q + j*q + 0] = 2 * ((er[i*x + j] - 1) * Ex[i*x + j]) - ex[i*x*q + j*q + 0];
				exb[i*x*q + j*q + 1] = 2 * W * (Ex[i*x + j] - V[1*3 + 1] * Hz[i*x + j] + V[1*3 + 2] * Hy[i*x + j]) - ex[i*x*q + j*q + 1];
				exb[i*x*q + j*q + 2] = 2 * W * (Ex[i*x + j] - V[2*3 + 1] * Hz[i*x + j] + V[2*3 + 2] * Hy[i*x + j]) - ex[i*x*q + j*q + 2];
				exb[i*x*q + j*q + 3] = 2 * W * (Ex[i*x + j] - V[3*3 + 1] * Hz[i*x + j] + V[3*3 + 2] * Hy[i*x + j]) - ex[i*x*q + j*q + 3];
				exb[i*x*q + j*q + 4] = 2 * W * (Ex[i*x + j] - V[4*3 + 1] * Hz[i*x + j] + V[4*3 + 2] * Hy[i*x + j]) - ex[i*x*q + j*q + 4];
				exb[i*x*q + j*q + 5] = 2 * W * (Ex[i*x + j] - V[5*3 + 1] * Hz[i*x + j] + V[5*3 + 2] * Hy[i*x + j]) - ex[i*x*q + j*q + 5];
				exb[i*x*q + j*q + 6] = 2 * W * (Ex[i*x + j] - V[6*3 + 1] * Hz[i*x + j] + V[6*3 + 2] * Hy[i*x + j]) - ex[i*x*q + j*q + 6];
				
				
				eyb[i*x*q + j*q + 0] = 2 * ((er[i*x + j] - 1) * Ey[i*x + j]) - ey[i*x*q + j*q + 0];
				eyb[i*x*q + j*q + 1] = 2 * W * (Ey[i*x + j] - V[1*3 + 2] * Hx[i*x + j] + V[1*3 + 0] * Hz[i*x + j]) - ey[i*x*q + j*q + 1];
				eyb[i*x*q + j*q + 2] = 2 * W * (Ey[i*x + j] - V[2*3 + 2] * Hx[i*x + j] + V[2*3 + 0] * Hz[i*x + j]) - ey[i*x*q + j*q + 2];
				eyb[i*x*q + j*q + 3] = 2 * W * (Ey[i*x + j] - V[3*3 + 2] * Hx[i*x + j] + V[3*3 + 0] * Hz[i*x + j]) - ey[i*x*q + j*q + 3];
				eyb[i*x*q + j*q + 4] = 2 * W * (Ey[i*x + j] - V[4*3 + 2] * Hx[i*x + j] + V[4*3 + 0] * Hz[i*x + j]) - ey[i*x*q + j*q + 4];
				eyb[i*x*q + j*q + 5] = 2 * W * (Ey[i*x + j] - V[5*3 + 2] * Hx[i*x + j] + V[5*3 + 0] * Hz[i*x + j]) - ey[i*x*q + j*q + 5];
				eyb[i*x*q + j*q + 6] = 2 * W * (Ey[i*x + j] - V[6*3 + 2] * Hx[i*x + j] + V[6*3 + 0] * Hz[i*x + j]) - ey[i*x*q + j*q + 6];


				ezb[i*x*q + j*q + 0] = 2 * ((er[i*x + j] - 1) * Ez[i*x + j]) - ez[i*x*q + j*q + 0];
				ezb[i*x*q + j*q + 1] = 2 * W * (Ez[i*x + j] - V[1*3 + 0] * Hy[i*x + j] + V[1*3 + 1] * Hx[i*x + j]) - ez[i*x*q + j*q + 1];
				ezb[i*x*q + j*q + 2] = 2 * W * (Ez[i*x + j] - V[2*3 + 0] * Hy[i*x + j] + V[2*3 + 1] * Hx[i*x + j]) - ez[i*x*q + j*q + 2];
				ezb[i*x*q + j*q + 3] = 2 * W * (Ez[i*x + j] - V[3*3 + 0] * Hy[i*x + j] + V[3*3 + 1] * Hx[i*x + j]) - ez[i*x*q + j*q + 3];
				ezb[i*x*q + j*q + 4] = 2 * W * (Ez[i*x + j] - V[4*3 + 0] * Hy[i*x + j] + V[4*3 + 1] * Hx[i*x + j]) - ez[i*x*q + j*q + 4];
				ezb[i*x*q + j*q + 5] = 2 * W * (Ez[i*x + j] - V[5*3 + 0] * Hy[i*x + j] + V[5*3 + 1] * Hx[i*x + j]) - ez[i*x*q + j*q + 5];
				ezb[i*x*q + j*q + 6] = 2 * W * (Ez[i*x + j] - V[6*3 + 0] * Hy[i*x + j] + V[6*3 + 1] * Hx[i*x + j]) - ez[i*x*q + j*q + 6];



				hxb[i*x*q + j*q + 0] = 2 * ((mur[i*x + j] - 1) * Hx[i*x + j]) - hx[i*x*q + j*q + 0];
				hxb[i*x*q + j*q + 1] = 2 * W * (Hx[i*x + j] - V[1*3 + 2] * Ey[i*x + j] + V[1*3 + 1] * Ez[i*x + j]) - hx[i*x*q + j*q + 1];
				hxb[i*x*q + j*q + 2] = 2 * W * (Hx[i*x + j] - V[2*3 + 2] * Ey[i*x + j] + V[2*3 + 1] * Ez[i*x + j]) - hx[i*x*q + j*q + 2];
				hxb[i*x*q + j*q + 3] = 2 * W * (Hx[i*x + j] - V[3*3 + 2] * Ey[i*x + j] + V[3*3 + 1] * Ez[i*x + j]) - hx[i*x*q + j*q + 3];
				hxb[i*x*q + j*q + 4] = 2 * W * (Hx[i*x + j] - V[4*3 + 2] * Ey[i*x + j] + V[4*3 + 1] * Ez[i*x + j]) - hx[i*x*q + j*q + 4];
				hxb[i*x*q + j*q + 5] = 2 * W * (Hx[i*x + j] - V[5*3 + 2] * Ey[i*x + j] + V[5*3 + 1] * Ez[i*x + j]) - hx[i*x*q + j*q + 5];
				hxb[i*x*q + j*q + 6] = 2 * W * (Hx[i*x + j] - V[6*3 + 2] * Ey[i*x + j] + V[6*3 + 1] * Ez[i*x + j]) - hx[i*x*q + j*q + 6];
				
				
				
				hyb[i*x*q + j*q + 0] = 2 * ((mur[i*x + j] - 1) * Hy[i*x + j]) - hy[i*x*q + j*q + 0];
				hyb[i*x*q + j*q + 1] = 2 * W * (Hy[i*x + j] - V[1*3 + 0] * Ez[i*x + j] + V[1*3 + 2] * Ex[i*x + j]) - hy[i*x*q + j*q + 1];
				hyb[i*x*q + j*q + 2] = 2 * W * (Hy[i*x + j] - V[2*3 + 0] * Ez[i*x + j] + V[2*3 + 2] * Ex[i*x + j]) - hy[i*x*q + j*q + 2];
				hyb[i*x*q + j*q + 3] = 2 * W * (Hy[i*x + j] - V[3*3 + 0] * Ez[i*x + j] + V[3*3 + 2] * Ex[i*x + j]) - hy[i*x*q + j*q + 3];
				hyb[i*x*q + j*q + 4] = 2 * W * (Hy[i*x + j] - V[4*3 + 0] * Ez[i*x + j] + V[4*3 + 2] * Ex[i*x + j]) - hy[i*x*q + j*q + 4];
				hyb[i*x*q + j*q + 5] = 2 * W * (Hy[i*x + j] - V[5*3 + 0] * Ez[i*x + j] + V[5*3 + 2] * Ex[i*x + j]) - hy[i*x*q + j*q + 5];
				hyb[i*x*q + j*q + 6] = 2 * W * (Hy[i*x + j] - V[6*3 + 0] * Ez[i*x + j] + V[6*3 + 2] * Ex[i*x + j]) - hy[i*x*q + j*q + 6];


				hzb[i*x*q + j*q + 0] = 2 * ((mur[i*x + j] - 1) * Hz[i*x + j]) - hz[i*x*q + j*q + 0];
				hzb[i*x*q + j*q + 1] = 2 * W * (Hz[i*x + j] - V[1*3 + 1] * Ex[i*x + j] + V[1*3 + 0] * Ey[i*x + j]) - hz[i*x*q + j*q + 1];
				hzb[i*x*q + j*q + 2] = 2 * W * (Hz[i*x + j] - V[2*3 + 1] * Ex[i*x + j] + V[2*3 + 0] * Ey[i*x + j]) - hz[i*x*q + j*q + 2];
				hzb[i*x*q + j*q + 3] = 2 * W * (Hz[i*x + j] - V[3*3 + 1] * Ex[i*x + j] + V[3*3 + 0] * Ey[i*x + j]) - hz[i*x*q + j*q + 3];
				hzb[i*x*q + j*q + 4] = 2 * W * (Hz[i*x + j] - V[4*3 + 1] * Ex[i*x + j] + V[4*3 + 0] * Ey[i*x + j]) - hz[i*x*q + j*q + 4];
				hzb[i*x*q + j*q + 5] = 2 * W * (Hz[i*x + j] - V[5*3 + 1] * Ex[i*x + j] + V[5*3 + 0] * Ey[i*x + j]) - hz[i*x*q + j*q + 5];
				hzb[i*x*q + j*q + 6] = 2 * W * (Hz[i*x + j] - V[6*3 + 1] * Ex[i*x + j] + V[6*3 + 0] * Ey[i*x + j]) - hz[i*x*q + j*q + 6];
			}
		}
	}
}




/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/










/*---------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------COLLISION + STREAMING IF FIELDS ARE NOT FORCED----------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/


void collNotForcingNode(float* ex, float* ey, float* ez, float* hx, float* hy, float* hz, float* exb, float* eyb, float* ezb, float* hxb, float* hyb, float* hzb, float* Ex, float* Ey, float* Ez, float* Hx, float* Hy, float* Hz, float* er, float* mur, int y, int x, int q) {

	int i, j;
	
	
		
	for (i = 0; i < y; i++) {
		omp_set_num_threads(N);
		#pragma omp parallel for
		for (j = 0; j < x; j++) {
				
			exb[i*x*q + j*q + 0] = 2 * ((er[i*x + j] - 1) * Ex[i*x + j]) - ex[i*x*q + j*q + 0];
			exb[i*x*q + j*q + 1] = 2 * W * (Ex[i*x + j] - V[1*3 + 1] * Hz[i*x + j] + V[1*3 + 2] * Hy[i*x + j]) - ex[i*x*q + j*q + 1];
			exb[i*x*q + j*q + 2] = 2 * W * (Ex[i*x + j] - V[2*3 + 1] * Hz[i*x + j] + V[2*3 + 2] * Hy[i*x + j]) - ex[i*x*q + j*q + 2];
			exb[i*x*q + j*q + 3] = 2 * W * (Ex[i*x + j] - V[3*3 + 1] * Hz[i*x + j] + V[3*3 + 2] * Hy[i*x + j]) - ex[i*x*q + j*q + 3];
			exb[i*x*q + j*q + 4] = 2 * W * (Ex[i*x + j] - V[4*3 + 1] * Hz[i*x + j] + V[4*3 + 2] * Hy[i*x + j]) - ex[i*x*q + j*q + 4];
			exb[i*x*q + j*q + 5] = 2 * W * (Ex[i*x + j] - V[5*3 + 1] * Hz[i*x + j] + V[5*3 + 2] * Hy[i*x + j]) - ex[i*x*q + j*q + 5];
			exb[i*x*q + j*q + 6] = 2 * W * (Ex[i*x + j] - V[6*3 + 1] * Hz[i*x + j] + V[6*3 + 2] * Hy[i*x + j]) - ex[i*x*q + j*q + 6];
			
			
			eyb[i*x*q + j*q + 0] = 2 * ((er[i*x + j] - 1) * Ey[i*x + j]) - ey[i*x*q + j*q + 0];
			eyb[i*x*q + j*q + 1] = 2 * W * (Ey[i*x + j] - V[1*3 + 2] * Hx[i*x + j] + V[1*3 + 0] * Hz[i*x + j]) - ey[i*x*q + j*q + 1];
			eyb[i*x*q + j*q + 2] = 2 * W * (Ey[i*x + j] - V[2*3 + 2] * Hx[i*x + j] + V[2*3 + 0] * Hz[i*x + j]) - ey[i*x*q + j*q + 2];
			eyb[i*x*q + j*q + 3] = 2 * W * (Ey[i*x + j] - V[3*3 + 2] * Hx[i*x + j] + V[3*3 + 0] * Hz[i*x + j]) - ey[i*x*q + j*q + 3];
			eyb[i*x*q + j*q + 4] = 2 * W * (Ey[i*x + j] - V[4*3 + 2] * Hx[i*x + j] + V[4*3 + 0] * Hz[i*x + j]) - ey[i*x*q + j*q + 4];
			eyb[i*x*q + j*q + 5] = 2 * W * (Ey[i*x + j] - V[5*3 + 2] * Hx[i*x + j] + V[5*3 + 0] * Hz[i*x + j]) - ey[i*x*q + j*q + 5];
			eyb[i*x*q + j*q + 6] = 2 * W * (Ey[i*x + j] - V[6*3 + 2] * Hx[i*x + j] + V[6*3 + 0] * Hz[i*x + j]) - ey[i*x*q + j*q + 6];


			ezb[i*x*q + j*q + 0] = 2 * ((er[i*x + j] - 1) * Ez[i*x + j]) - ez[i*x*q + j*q + 0];
			ezb[i*x*q + j*q + 1] = 2 * W * (Ez[i*x + j] - V[1*3 + 0] * Hy[i*x + j] + V[1*3 + 1] * Hx[i*x + j]) - ez[i*x*q + j*q + 1];
			ezb[i*x*q + j*q + 2] = 2 * W * (Ez[i*x + j] - V[2*3 + 0] * Hy[i*x + j] + V[2*3 + 1] * Hx[i*x + j]) - ez[i*x*q + j*q + 2];
			ezb[i*x*q + j*q + 3] = 2 * W * (Ez[i*x + j] - V[3*3 + 0] * Hy[i*x + j] + V[3*3 + 1] * Hx[i*x + j]) - ez[i*x*q + j*q + 3];
			ezb[i*x*q + j*q + 4] = 2 * W * (Ez[i*x + j] - V[4*3 + 0] * Hy[i*x + j] + V[4*3 + 1] * Hx[i*x + j]) - ez[i*x*q + j*q + 4];
			ezb[i*x*q + j*q + 5] = 2 * W * (Ez[i*x + j] - V[5*3 + 0] * Hy[i*x + j] + V[5*3 + 1] * Hx[i*x + j]) - ez[i*x*q + j*q + 5];
			ezb[i*x*q + j*q + 6] = 2 * W * (Ez[i*x + j] - V[6*3 + 0] * Hy[i*x + j] + V[6*3 + 1] * Hx[i*x + j]) - ez[i*x*q + j*q + 6];



			hxb[i*x*q + j*q + 0] = 2 * ((mur[i*x + j] - 1) * Hx[i*x + j]) - hx[i*x*q + j*q + 0];
			hxb[i*x*q + j*q + 1] = 2 * W * (Hx[i*x + j] - V[1*3 + 2] * Ey[i*x + j] + V[1*3 + 1] * Ez[i*x + j]) - hx[i*x*q + j*q + 1];
			hxb[i*x*q + j*q + 2] = 2 * W * (Hx[i*x + j] - V[2*3 + 2] * Ey[i*x + j] + V[2*3 + 1] * Ez[i*x + j]) - hx[i*x*q + j*q + 2];
			hxb[i*x*q + j*q + 3] = 2 * W * (Hx[i*x + j] - V[3*3 + 2] * Ey[i*x + j] + V[3*3 + 1] * Ez[i*x + j]) - hx[i*x*q + j*q + 3];
			hxb[i*x*q + j*q + 4] = 2 * W * (Hx[i*x + j] - V[4*3 + 2] * Ey[i*x + j] + V[4*3 + 1] * Ez[i*x + j]) - hx[i*x*q + j*q + 4];
			hxb[i*x*q + j*q + 5] = 2 * W * (Hx[i*x + j] - V[5*3 + 2] * Ey[i*x + j] + V[5*3 + 1] * Ez[i*x + j]) - hx[i*x*q + j*q + 5];
			hxb[i*x*q + j*q + 6] = 2 * W * (Hx[i*x + j] - V[6*3 + 2] * Ey[i*x + j] + V[6*3 + 1] * Ez[i*x + j]) - hx[i*x*q + j*q + 6];
			
			
			
			hyb[i*x*q + j*q + 0] = 2 * ((mur[i*x + j] - 1) * Hy[i*x + j]) - hy[i*x*q + j*q + 0];
			hyb[i*x*q + j*q + 1] = 2 * W * (Hy[i*x + j] - V[1*3 + 0] * Ez[i*x + j] + V[1*3 + 2] * Ex[i*x + j]) - hy[i*x*q + j*q + 1];
			hyb[i*x*q + j*q + 2] = 2 * W * (Hy[i*x + j] - V[2*3 + 0] * Ez[i*x + j] + V[2*3 + 2] * Ex[i*x + j]) - hy[i*x*q + j*q + 2];
			hyb[i*x*q + j*q + 3] = 2 * W * (Hy[i*x + j] - V[3*3 + 0] * Ez[i*x + j] + V[3*3 + 2] * Ex[i*x + j]) - hy[i*x*q + j*q + 3];
			hyb[i*x*q + j*q + 4] = 2 * W * (Hy[i*x + j] - V[4*3 + 0] * Ez[i*x + j] + V[4*3 + 2] * Ex[i*x + j]) - hy[i*x*q + j*q + 4];
			hyb[i*x*q + j*q + 5] = 2 * W * (Hy[i*x + j] - V[5*3 + 0] * Ez[i*x + j] + V[5*3 + 2] * Ex[i*x + j]) - hy[i*x*q + j*q + 5];
			hyb[i*x*q + j*q + 6] = 2 * W * (Hy[i*x + j] - V[6*3 + 0] * Ez[i*x + j] + V[6*3 + 2] * Ex[i*x + j]) - hy[i*x*q + j*q + 6];


			hzb[i*x*q + j*q + 0] = 2 * ((mur[i*x + j] - 1) * Hz[i*x + j]) - hz[i*x*q + j*q + 0];
			hzb[i*x*q + j*q + 1] = 2 * W * (Hz[i*x + j] - V[1*3 + 1] * Ex[i*x + j] + V[1*3 + 0] * Ey[i*x + j]) - hz[i*x*q + j*q + 1];
			hzb[i*x*q + j*q + 2] = 2 * W * (Hz[i*x + j] - V[2*3 + 1] * Ex[i*x + j] + V[2*3 + 0] * Ey[i*x + j]) - hz[i*x*q + j*q + 2];
			hzb[i*x*q + j*q + 3] = 2 * W * (Hz[i*x + j] - V[3*3 + 1] * Ex[i*x + j] + V[3*3 + 0] * Ey[i*x + j]) - hz[i*x*q + j*q + 3];
			hzb[i*x*q + j*q + 4] = 2 * W * (Hz[i*x + j] - V[4*3 + 1] * Ex[i*x + j] + V[4*3 + 0] * Ey[i*x + j]) - hz[i*x*q + j*q + 4];
			hzb[i*x*q + j*q + 5] = 2 * W * (Hz[i*x + j] - V[5*3 + 1] * Ex[i*x + j] + V[5*3 + 0] * Ey[i*x + j]) - hz[i*x*q + j*q + 5];
			hzb[i*x*q + j*q + 6] = 2 * W * (Hz[i*x + j] - V[6*3 + 1] * Ex[i*x + j] + V[6*3 + 0] * Ey[i*x + j]) - hz[i*x*q + j*q + 6];
			
		}
	}
}


/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/










/*---------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------COLLISION + STREAMING IF FIELDS ARE NOT FORCED----------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/


void streaming(float* ex, float* ey, float* ez, float* hx, float* hy, float* hz, float* exb, float* eyb, float* ezb, float* hxb, float* hyb, float* hzb, int y, int x, int q) {

	int i, j;
	
	
		
	for (i = 0; i < y; i++) {
		omp_set_num_threads(N);
		#pragma omp parallel for
		for (j = 0; j < x; j++) {
		
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
				
			ex[i*x*q + j*q + 0] = exb[i*x*q + j*q  + 0];
			ex[i*x*q + j*q + 1] = exb[i*x*q + j1*q + 1];
			ex[i*x*q + j*q + 2] = exb[i*x*q + j2*q + 2];
			ex[i*x*q + j*q + 3] = exb[i3*x*q + j*q + 3];
			ex[i*x*q + j*q + 4] = exb[i4*x*q + j*q + 4];
			ex[i*x*q + j*q + 5] = exb[i*x*q + j*q  + 5];
			ex[i*x*q + j*q + 6] = exb[i*x*q + j*q  + 6];


			ey[i*x*q + j*q + 0] = eyb[i*x*q + j*q  + 0];
			ey[i*x*q + j*q + 1] = eyb[i*x*q + j1*q + 1];
			ey[i*x*q + j*q + 2] = eyb[i*x*q + j2*q + 2];
			ey[i*x*q + j*q + 3] = eyb[i3*x*q + j*q + 3];
			ey[i*x*q + j*q + 4] = eyb[i4*x*q + j*q + 4];
			ey[i*x*q + j*q + 5] = eyb[i*x*q + j*q  + 5];
			ey[i*x*q + j*q + 6] = eyb[i*x*q + j*q  + 6];
			
			
			ez[i*x*q + j*q + 0] = ezb[i*x*q + j*q  + 0];
			ez[i*x*q + j*q + 1] = ezb[i*x*q + j1*q + 1];
			ez[i*x*q + j*q + 2] = ezb[i*x*q + j2*q + 2];
			ez[i*x*q + j*q + 3] = ezb[i3*x*q + j*q + 3];
			ez[i*x*q + j*q + 4] = ezb[i4*x*q + j*q + 4];
			ez[i*x*q + j*q + 5] = ezb[i*x*q + j*q  + 5];
			ez[i*x*q + j*q + 6] = ezb[i*x*q + j*q  + 6];
			



			hx[i*x*q + j*q + 0] = hxb[i*x*q + j*q  + 0];
			hx[i*x*q + j*q + 1] = hxb[i*x*q + j1*q + 1];
			hx[i*x*q + j*q + 2] = hxb[i*x*q + j2*q + 2];
			hx[i*x*q + j*q + 3] = hxb[i3*x*q + j*q + 3];
			hx[i*x*q + j*q + 4] = hxb[i4*x*q + j*q + 4];
			hx[i*x*q + j*q + 5] = hxb[i*x*q + j*q  + 5];
			hx[i*x*q + j*q + 6] = hxb[i*x*q + j*q  + 6];


			hy[i*x*q + j*q + 0] = hyb[i*x*q + j*q  + 0];
			hy[i*x*q + j*q + 1] = hyb[i*x*q + j1*q + 1];
			hy[i*x*q + j*q + 2] = hyb[i*x*q + j2*q + 2];
			hy[i*x*q + j*q + 3] = hyb[i3*x*q + j*q + 3];
			hy[i*x*q + j*q + 4] = hyb[i4*x*q + j*q + 4];
			hy[i*x*q + j*q + 5] = hyb[i*x*q + j*q  + 5];
			hy[i*x*q + j*q + 6] = hyb[i*x*q + j*q  + 6];
			
			
			hz[i*x*q + j*q + 0] = hzb[i*x*q + j*q  + 0];
			hz[i*x*q + j*q + 1] = hzb[i*x*q + j1*q + 1];
			hz[i*x*q + j*q + 2] = hzb[i*x*q + j2*q + 2];
			hz[i*x*q + j*q + 3] = hzb[i3*x*q + j*q + 3];
			hz[i*x*q + j*q + 4] = hzb[i4*x*q + j*q + 4];
			hz[i*x*q + j*q + 5] = hzb[i*x*q + j*q  + 5];
			hz[i*x*q + j*q + 6] = hzb[i*x*q + j*q  + 6];
			

		}
	}
}


/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/







