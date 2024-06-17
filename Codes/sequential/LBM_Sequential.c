#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define W 1.0/6.0



int V[21] = {0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1};


/*---------------------------------------------------------------------------------------------*/
// calculates the macroscopic fields from the distribution functions


void macroField(float* dis_func, float* mat_prop, float* field, int z, int y, int x, int q) {
	
	int i, j, k, l;
	
	for(i = 0; i < z; i++) {
		for(j = 0; j < y; j++) {
			for(k = 0; k < x; k++) {
				for(l = 0; l < q; l++) {
					field[i*y*x + j*x + k] += dis_func[i*y*x*q + j*x*q + k*q + l] / mat_prop[i*y*x + j*x + k];
				}
			}
		}
	}
}
/*--------------------------------------------------------------------------------------------*/






/*----------------------------------------------------------------------------------------------------*/
// inintialize the macroscopic field


void initializeField(float* Ex, float* Ey, float* Ez, float* Hx, float* Hy, float* Hz, int z, int y, int x) {
	
	for(int i = 0; i < z; i++){
		for(int j = 0; j < y; j++){
			for(int k = 0; k < x; k++){				
				Ex[i*y*x + j*x + k] = 0;
				Ey[i*y*x + j*x + k] = 0;
				Ez[i*y*x + j*x + k] = 0;
				
				Hx[i*y*x + j*x + k] = 0;
				Hy[i*y*x + j*x + k] = 0;
				Hz[i*y*x + j*x + k] = 0;
			}
		}
	}
}

/*----------------------------------------------------------------------------------------------------*/








/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------COLLISION + STREAMING IF FIELDS ARE FORCED------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/


void collStreamForcingNode(float* ex, float* ey, float* ez, float* hx, float* hy, float* hz, float* Ex, float* Ey, float* Ez, float* Hx, float* Hy, float* Hz, float* er, float* mur, int z, int y, int x, int q, int xmin, int xmax, int ymin, int ymax) {

	int i=0, j=0, k=0;
	
	

	for (i = 0; i < z; i++) {
		for (j = 0; j < y; j++) {
			for (k = 0; k < x; k++) {
			
				int k1 = x-1-k, k11 = k1-1;
				int k21 = k+1;
				int j3 = y-1-j, j31 = j3-1;
				int j41 = j+1;

				if (k1 == 0)
					k11 = k1;
				if (k == x-1)
					k21 = k;
				if (j3 == 0)
					j31 = j3;
				if (j == y-1)
					j41 = j;
				
				if (j >= ymin && j <= ymax && k >= xmin && k <= xmax) { // no collision, only streaming at the forcing nodes

					ex[i*y*x*q + j*x*q + k*q  + 0] = ((er[i*y*x + j*x + k] - 1) * Ex[i*y*x + j*x + k]);	
					ex[i*y*x*q + j*x*q + k1*q + 1] = W * (Ex[i*y*x + j*x + k11] - V[1*3 + 1] * Hz[i*y*x + j*x + k11] + V[1*3 + 2] * Hy[i*y*x + j*x + k11]);
					ex[i*y*x*q + j*x*q + k*q  + 2] = W * (Ex[i*y*x + j*x + k21] - V[2*3 + 1] * Hz[i*y*x + j*x + k21] + V[2*3 + 2] * Hy[i*y*x + j*x + k21]);
					ex[i*y*x*q + j3*x*q + k*q + 3] = W * (Ex[i*y*x + j31*x + k] - V[3*3 + 1] * Hz[i*y*x + j31*x + k] + V[3*3 + 2] * Hy[i*y*x + j31*x + k]);
					ex[i*y*x*q + j*x*q + k*q  + 4] = W * (Ex[i*y*x + j41*x + k] - V[4*3 + 1] * Hz[i*y*x + j41*x + k] + V[4*3 + 2] * Hy[i*y*x + j41*x + k]);
					ex[i*y*x*q + j*x*q + k*q  + 5] = W * (Ex[i*y*x + j*x   + k] - V[5*3 + 1] * Hz[i*y*x + j*x   + k] + V[5*3 + 2] * Hy[i*y*x + j*x   + k]);
					ex[i*y*x*q + j*x*q + k*q  + 6] = W * (Ex[i*y*x + j*x   + k] - V[6*3 + 1] * Hz[i*y*x + j*x   + k] + V[6*3 + 2] * Hy[i*y*x + j*x   + k]);
					
					
					ey[i*y*x*q + j*x*q + k*q  + 0] = ((er[i*y*x + j*x + k] - 1) * Ey[i*y*x + j*x + k]);
					ey[i*y*x*q + j*x*q + k1*q + 1] = W * (Ey[i*y*x + j*x + k11] - V[1*3 + 2] * Hx[i*y*x + j*x + k11] + V[1*3 + 0] * Hz[i*y*x + j*x + k11]);
					ey[i*y*x*q + j*x*q + k*q  + 2] = W * (Ey[i*y*x + j*x + k21] - V[2*3 + 2] * Hx[i*y*x + j*x + k21] + V[2*3 + 0] * Hz[i*y*x + j*x + k21]);
					ey[i*y*x*q + j3*x*q + k*q + 3] = W * (Ey[i*y*x + j31*x + k] - V[3*3 + 2] * Hx[i*y*x + j31*x + k] + V[3*3 + 0] * Hz[i*y*x + j31*x + k]);
					ey[i*y*x*q + j*x*q + k*q  + 4] = W * (Ey[i*y*x + j41*x + k] - V[4*3 + 2] * Hx[i*y*x + j41*x + k] + V[4*3 + 0] * Hz[i*y*x + j41*x + k]);
					ey[i*y*x*q + j*x*q + k*q  + 5] = W * (Ey[i*y*x + j*x   + k] - V[5*3 + 2] * Hx[i*y*x + j*x   + k] + V[5*3 + 0] * Hz[i*y*x + j*x   + k]);
					ey[i*y*x*q + j*x*q + k*q  + 6] = W * (Ey[i*y*x + j*x   + k] - V[6*3 + 2] * Hx[i*y*x + j*x   + k] + V[6*3 + 0] * Hz[i*y*x + j*x   + k]);


					ez[i*y*x*q + j*x*q + k*q  + 0] = ((er[i*y*x + j*x + k] - 1) * Ez[i*y*x + j*x + k]);
					ez[i*y*x*q + j*x*q + k1*q + 1] = W * (Ez[i*y*x + j*x + k11] - V[1*3 + 0] * Hy[i*y*x + j*x + k11] + V[1*3 + 1] * Hx[i*y*x + j*x + k11]);
					ez[i*y*x*q + j*x*q + k*q  + 2] = W * (Ez[i*y*x + j*x + k21] - V[2*3 + 0] * Hy[i*y*x + j*x + k21] + V[2*3 + 1] * Hx[i*y*x + j*x + k21]);
					ez[i*y*x*q + j3*x*q + k*q + 3] = W * (Ez[i*y*x + j31*x + k] - V[3*3 + 0] * Hy[i*y*x + j31*x + k] + V[3*3 + 1] * Hx[i*y*x + j31*x + k]);
					ez[i*y*x*q + j*x*q + k*q  + 4] = W * (Ez[i*y*x + j41*x + k] - V[4*3 + 0] * Hy[i*y*x + j41*x + k] + V[4*3 + 1] * Hx[i*y*x + j41*x + k]);
					ez[i*y*x*q + j*x*q + k*q  + 5] = W * (Ez[i*y*x + j*x + k]   - V[5*3 + 0] * Hy[i*y*x + j*x   + k] + V[5*3 + 1] * Hx[i*y*x + j*x   + k]);
					ez[i*y*x*q + j*x*q + k*q  + 6] = W * (Ez[i*y*x + j*x + k]   - V[6*3 + 0] * Hy[i*y*x + j*x   + k] + V[6*3 + 1] * Hx[i*y*x + j*x   + k]);



					hx[i*y*x*q + j*x*q + k*q  + 0] = ((mur[i*y*x + j*x + k] - 1) * Hx[i*y*x + j*x + k]);
					hx[i*y*x*q + j*x*q + k1*q + 1] = W * (Hx[i*y*x + j*x + k11] - V[1*3 + 2] * Ey[i*y*x + j*x + k11] + V[1*3 + 1] * Ez[i*y*x + j*x + k11]);
					hx[i*y*x*q + j*x*q + k*q  + 2] = W * (Hx[i*y*x + j*x + k21] - V[2*3 + 2] * Ey[i*y*x + j*x + k21] + V[2*3 + 1] * Ez[i*y*x + j*x + k21]);
					hx[i*y*x*q + j3*x*q + k*q + 3] = W * (Hx[i*y*x + j31*x + k] - V[3*3 + 2] * Ey[i*y*x + j31*x + k] + V[3*3 + 1] * Ez[i*y*x + j31*x + k]);
					hx[i*y*x*q + j*x*q + k*q  + 4] = W * (Hx[i*y*x + j41*x + k] - V[4*3 + 2] * Ey[i*y*x + j41*x + k] + V[4*3 + 1] * Ez[i*y*x + j41*x + k]);
					hx[i*y*x*q + j*x*q + k*q  + 5] = W * (Hx[i*y*x + j*x   + k] - V[5*3 + 2] * Ey[i*y*x + j*x   + k] + V[5*3 + 1] * Ez[i*y*x + j*x   + k]);
					hx[i*y*x*q + j*x*q + k*q  + 6] = W * (Hx[i*y*x + j*x   + k] - V[6*3 + 2] * Ey[i*y*x + j*x   + k] + V[6*3 + 1] * Ez[i*y*x + j*x   + k]);
					
					
				
					hy[i*y*x*q + j*x*q + k*q  + 0] = ((mur[i*y*x + j*x + k] - 1) * Hy[i*y*x + j*x + k]);
					hy[i*y*x*q + j*x*q + k1*q + 1] = W * (Hy[i*y*x + j*x + k11] - V[1*3 + 0] * Ez[i*y*x + j*x + k11] + V[1*3 + 2] * Ex[i*y*x + j*x + k11]);
					hy[i*y*x*q + j*x*q + k*q  + 2] = W * (Hy[i*y*x + j*x + k21] - V[2*3 + 0] * Ez[i*y*x + j*x + k21] + V[2*3 + 2] * Ex[i*y*x + j*x + k21]);
					hy[i*y*x*q + j3*x*q + k*q + 3] = W * (Hy[i*y*x + j31*x + k] - V[3*3 + 0] * Ez[i*y*x + j31*x + k] + V[3*3 + 2] * Ex[i*y*x + j31*x + k]);
					hy[i*y*x*q + j*x*q + k*q  + 4] = W * (Hy[i*y*x + j41*x + k] - V[4*3 + 0] * Ez[i*y*x + j41*x + k] + V[4*3 + 2] * Ex[i*y*x + j41*x + k]);
					hy[i*y*x*q + j*x*q + k*q  + 5] = W * (Hy[i*y*x + j*x   + k] - V[5*3 + 0] * Ez[i*y*x + j*x   + k] + V[5*3 + 2] * Ex[i*y*x + j*x   + k]);
					hy[i*y*x*q + j*x*q + k*q  + 6] = W * (Hy[i*y*x + j*x   + k] - V[6*3 + 0] * Ez[i*y*x + j*x   + k] + V[6*3 + 2] * Ex[i*y*x + j*x   + k]);


					hz[i*y*x*q + j*x*q + k*q  + 0] = ((mur[i*y*x + j*x + k] - 1) * Hz[i*y*x + j*x + k]);
					hz[i*y*x*q + j*x*q + k1*q + 1] = W * (Hz[i*y*x + j*x + k11] - V[1*3 + 1] * Ex[i*y*x + j*x + k11] + V[1*3 + 0] * Ey[i*y*x + j*x + k11]);
					hz[i*y*x*q + j*x*q + k*q  + 2] = W * (Hz[i*y*x + j*x + k21] - V[2*3 + 1] * Ex[i*y*x + j*x + k21] + V[2*3 + 0] * Ey[i*y*x + j*x + k21]);
					hz[i*y*x*q + j3*x*q + k*q + 3] = W * (Hz[i*y*x + j31*x + k] - V[3*3 + 1] * Ex[i*y*x + j31*x + k] + V[3*3 + 0] * Ey[i*y*x + j31*x + k]);
					hz[i*y*x*q + j*x*q + k*q  + 4] = W * (Hz[i*y*x + j41*x + k] - V[4*3 + 1] * Ex[i*y*x + j41*x + k] + V[4*3 + 0] * Ey[i*y*x + j41*x + k]);
					hz[i*y*x*q + j*x*q + k*q  + 5] = W * (Hz[i*y*x + j*x   + k] - V[5*3 + 1] * Ex[i*y*x + j*x   + k] + V[5*3 + 0] * Ey[i*y*x + j*x   + k]);
					hz[i*y*x*q + j*x*q + k*q  + 6] = W * (Hz[i*y*x + j*x   + k] - V[6*3 + 1] * Ex[i*y*x + j*x   + k] + V[6*3 + 0] * Ey[i*y*x + j*x   + k]);
				}
				
								
				
				else {
					
					ex[i*y*x*q + j*x*q + k*q  + 0] = 2 * ((er[i*y*x + j*x + k] - 1) * Ex[i*y*x + j*x + k]) - ex[i*y*x*q + j*x*q + k*q + 0];
					ex[i*y*x*q + j*x*q + k1*q + 1] = 2 * W * (Ex[i*y*x + j*x + k11] - V[1*3 + 1] * Hz[i*y*x + j*x + k11] + V[1*3 + 2] * Hy[i*y*x + j*x + k11]) - ex[i*y*x*q + j*x*q + k11*q + 1];
					ex[i*y*x*q + j*x*q + k*q  + 2] = 2 * W * (Ex[i*y*x + j*x + k21] - V[2*3 + 1] * Hz[i*y*x + j*x + k21] + V[2*3 + 2] * Hy[i*y*x + j*x + k21]) - ex[i*y*x*q + j*x*q + k21*q + 2];
					ex[i*y*x*q + j3*x*q + k*q + 3] = 2 * W * (Ex[i*y*x + j31*x + k] - V[3*3 + 1] * Hz[i*y*x + j31*x + k] + V[3*3 + 2] * Hy[i*y*x + j31*x + k]) - ex[i*y*x*q + j31*x*q + k*q + 3];
					ex[i*y*x*q + j*x*q + k*q  + 4] = 2 * W * (Ex[i*y*x + j41*x + k] - V[4*3 + 1] * Hz[i*y*x + j41*x + k] + V[4*3 + 2] * Hy[i*y*x + j41*x + k]) - ex[i*y*x*q + j41*x*q + k*q + 4];
					ex[i*y*x*q + j*x*q + k*q  + 5] = 2 * W * (Ex[i*y*x + j*x   + k] - V[5*3 + 1] * Hz[i*y*x + j*x   + k] + V[5*3 + 2] * Hy[i*y*x + j*x   + k]) - ex[i*y*x*q + j*x*q + k*q   + 5];
					ex[i*y*x*q + j*x*q + k*q  + 6] = 2 * W * (Ex[i*y*x + j*x   + k] - V[6*3 + 1] * Hz[i*y*x + j*x   + k] + V[6*3 + 2] * Hy[i*y*x + j*x   + k]) - ex[i*y*x*q + j*x*q + k*q   + 6];
					
					
					ey[i*y*x*q + j*x*q + k*q  + 0] = 2 * ((er[i*y*x + j*x + k] - 1) * Ey[i*y*x + j*x + k]) - ey[i*y*x*q + j*x*q + k*q + 0];
					ey[i*y*x*q + j*x*q + k1*q + 1] = 2 * W * (Ey[i*y*x + j*x + k11] - V[1*3 + 2] * Hx[i*y*x + j*x + k11] + V[1*3 + 0] * Hz[i*y*x + j*x + k11]) - ey[i*y*x*q + j*x*q + k11*q + 1];
					ey[i*y*x*q + j*x*q + k*q  + 2] = 2 * W * (Ey[i*y*x + j*x + k21] - V[2*3 + 2] * Hx[i*y*x + j*x + k21] + V[2*3 + 0] * Hz[i*y*x + j*x + k21]) - ey[i*y*x*q + j*x*q + k21*q + 2];
					ey[i*y*x*q + j3*x*q + k*q + 3] = 2 * W * (Ey[i*y*x + j31*x + k] - V[3*3 + 2] * Hx[i*y*x + j31*x + k] + V[3*3 + 0] * Hz[i*y*x + j31*x + k]) - ey[i*y*x*q + j31*x*q + k*q + 3];
					ey[i*y*x*q + j*x*q + k*q  + 4] = 2 * W * (Ey[i*y*x + j41*x + k] - V[4*3 + 2] * Hx[i*y*x + j41*x + k] + V[4*3 + 0] * Hz[i*y*x + j41*x + k]) - ey[i*y*x*q + j41*x*q + k*q + 4];
					ey[i*y*x*q + j*x*q + k*q  + 5] = 2 * W * (Ey[i*y*x + j*x   + k] - V[5*3 + 2] * Hx[i*y*x + j*x   + k] + V[5*3 + 0] * Hz[i*y*x + j*x   + k]) - ey[i*y*x*q + j*x*q + k*q   + 5];
					ey[i*y*x*q + j*x*q + k*q  + 6] = 2 * W * (Ey[i*y*x + j*x   + k] - V[6*3 + 2] * Hx[i*y*x + j*x   + k] + V[6*3 + 0] * Hz[i*y*x + j*x   + k]) - ey[i*y*x*q + j*x*q + k*q   + 6];


					ez[i*y*x*q + j*x*q + k*q  + 0] = 2 * ((er[i*y*x + j*x + k] - 1) * Ez[i*y*x + j*x + k]) - ez[i*y*x*q + j*x*q + k*q + 0];
					ez[i*y*x*q + j*x*q + k1*q + 1] = 2 * W * (Ez[i*y*x + j*x + k11] - V[1*3 + 0] * Hy[i*y*x + j*x + k11] + V[1*3 + 1] * Hx[i*y*x + j*x + k11]) - ez[i*y*x*q + j*x*q + k11*q + 1];
					ez[i*y*x*q + j*x*q + k*q  + 2] = 2 * W * (Ez[i*y*x + j*x + k21] - V[2*3 + 0] * Hy[i*y*x + j*x + k21] + V[2*3 + 1] * Hx[i*y*x + j*x + k21]) - ez[i*y*x*q + j*x*q + k21*q + 2];
					ez[i*y*x*q + j3*x*q + k*q + 3] = 2 * W * (Ez[i*y*x + j31*x + k] - V[3*3 + 0] * Hy[i*y*x + j31*x + k] + V[3*3 + 1] * Hx[i*y*x + j31*x + k]) - ez[i*y*x*q + j31*x*q + k*q + 3];
					ez[i*y*x*q + j*x*q + k*q  + 4] = 2 * W * (Ez[i*y*x + j41*x + k] - V[4*3 + 0] * Hy[i*y*x + j41*x + k] + V[4*3 + 1] * Hx[i*y*x + j41*x + k]) - ez[i*y*x*q + j41*x*q + k*q + 4];
					ez[i*y*x*q + j*x*q + k*q  + 5] = 2 * W * (Ez[i*y*x + j*x   + k] - V[5*3 + 0] * Hy[i*y*x + j*x   + k] + V[5*3 + 1] * Hx[i*y*x + j*x   + k]) - ez[i*y*x*q + j*x*q + k*q   + 5];
					ez[i*y*x*q + j*x*q + k*q  + 6] = 2 * W * (Ez[i*y*x + j*x   + k] - V[6*3 + 0] * Hy[i*y*x + j*x   + k] + V[6*3 + 1] * Hx[i*y*x + j*x   + k]) - ez[i*y*x*q + j*x*q + k*q   + 6];


					hx[i*y*x*q + j*x*q + k*q  + 0] = 2 * ((mur[i*y*x + j*x + k] - 1) * Hx[i*y*x + j*x + k]) - hx[i*y*x*q + j*x*q + k*q + 0];
					hx[i*y*x*q + j*x*q + k1*q + 1] = 2 * W * (Hx[i*y*x + j*x + k11] - V[1*3 + 2] * Ey[i*y*x + j*x + k11] + V[1*3 + 1] * Ez[i*y*x + j*x + k11]) - hx[i*y*x*q + j*x*q + k11*q + 1];
					hx[i*y*x*q + j*x*q + k*q  + 2] = 2 * W * (Hx[i*y*x + j*x + k21] - V[2*3 + 2] * Ey[i*y*x + j*x + k21] + V[2*3 + 1] * Ez[i*y*x + j*x + k21]) - hx[i*y*x*q + j*x*q + k21*q + 2];
					hx[i*y*x*q + j3*x*q + k*q + 3] = 2 * W * (Hx[i*y*x + j31*x + k] - V[3*3 + 2] * Ey[i*y*x + j31*x + k] + V[3*3 + 1] * Ez[i*y*x + j31*x + k]) - hx[i*y*x*q + j31*x*q + k*q + 3];
					hx[i*y*x*q + j*x*q + k*q  + 4] = 2 * W * (Hx[i*y*x + j41*x + k] - V[4*3 + 2] * Ey[i*y*x + j41*x + k] + V[4*3 + 1] * Ez[i*y*x + j41*x + k]) - hx[i*y*x*q + j41*x*q + k*q + 4];
					hx[i*y*x*q + j*x*q + k*q  + 5] = 2 * W * (Hx[i*y*x + j*x   + k] - V[5*3 + 2] * Ey[i*y*x + j*x   + k] + V[5*3 + 1] * Ez[i*y*x + j*x   + k]) - hx[i*y*x*q + j*x*q + k*q   + 5];
					hx[i*y*x*q + j*x*q + k*q  + 6] = 2 * W * (Hx[i*y*x + j*x   + k] - V[6*3 + 2] * Ey[i*y*x + j*x   + k] + V[6*3 + 1] * Ez[i*y*x + j*x   + k]) - hx[i*y*x*q + j*x*q + k*q   + 6];
					
					
					hy[i*y*x*q + j*x*q + k*q  + 0] = 2 * ((mur[i*y*x + j*x + k] - 1) * Hy[i*y*x + j*x + k]) - hy[i*y*x*q + j*x*q + k*q + 0];
					hy[i*y*x*q + j*x*q + k1*q + 1] = 2 * W * (Hy[i*y*x + j*x + k11] - V[1*3 + 0] * Ez[i*y*x + j*x + k11] + V[1*3 + 2] * Ex[i*y*x + j*x + k11]) - hy[i*y*x*q + j*x*q + k11*q + 1];
					hy[i*y*x*q + j*x*q + k*q  + 2] = 2 * W * (Hy[i*y*x + j*x + k21] - V[2*3 + 0] * Ez[i*y*x + j*x + k21] + V[2*3 + 2] * Ex[i*y*x + j*x + k21]) - hy[i*y*x*q + j*x*q + k21*q + 2];
					hy[i*y*x*q + j3*x*q + k*q + 3] = 2 * W * (Hy[i*y*x + j31*x + k] - V[3*3 + 0] * Ez[i*y*x + j31*x + k] + V[3*3 + 2] * Ex[i*y*x + j31*x + k]) - hy[i*y*x*q + j31*x*q + k*q + 3];
					hy[i*y*x*q + j*x*q + k*q  + 4] = 2 * W * (Hy[i*y*x + j41*x + k] - V[4*3 + 0] * Ez[i*y*x + j41*x + k] + V[4*3 + 2] * Ex[i*y*x + j41*x + k]) - hy[i*y*x*q + j41*x*q + k*q + 4];
					hy[i*y*x*q + j*x*q + k*q  + 5] = 2 * W * (Hy[i*y*x + j*x   + k] - V[5*3 + 0] * Ez[i*y*x + j*x   + k] + V[5*3 + 2] * Ex[i*y*x + j*x   + k]) - hy[i*y*x*q + j*x*q + k*q   + 5];
					hy[i*y*x*q + j*x*q + k*q  + 6] = 2 * W * (Hy[i*y*x + j*x   + k] - V[6*3 + 0] * Ez[i*y*x + j*x   + k] + V[6*3 + 2] * Ex[i*y*x + j*x   + k]) - hy[i*y*x*q + j*x*q + k*q   + 6];


					hz[i*y*x*q + j*x*q + k*q  + 0] = 2 * ((mur[i*y*x + j*x + k] - 1) * Hz[i*y*x + j*x + k]) - hz[i*y*x*q + j*x*q + k*q + 0];
					hz[i*y*x*q + j*x*q + k1*q + 1] = 2 * W * (Hz[i*y*x + j*x + k11] - V[1*3 + 1] * Ex[i*y*x + j*x + k11] + V[1*3 + 0] * Ey[i*y*x + j*x + k11]) - hz[i*y*x*q + j*x*q + k11*q + 1];
					hz[i*y*x*q + j*x*q + k*q  + 2] = 2 * W * (Hz[i*y*x + j*x + k21] - V[2*3 + 1] * Ex[i*y*x + j*x + k21] + V[2*3 + 0] * Ey[i*y*x + j*x + k21]) - hz[i*y*x*q + j*x*q + k21*q + 2];
					hz[i*y*x*q + j3*x*q + k*q + 3] = 2 * W * (Hz[i*y*x + j31*x + k] - V[3*3 + 1] * Ex[i*y*x + j31*x + k] + V[3*3 + 0] * Ey[i*y*x + j31*x + k]) - hz[i*y*x*q + j31*x*q + k*q + 3];
					hz[i*y*x*q + j*x*q + k*q  + 4] = 2 * W * (Hz[i*y*x + j41*x + k] - V[4*3 + 1] * Ex[i*y*x + j41*x + k] + V[4*3 + 0] * Ey[i*y*x + j41*x + k]) - hz[i*y*x*q + j41*x*q + k*q + 4];
					hz[i*y*x*q + j*x*q + k*q  + 5] = 2 * W * (Hz[i*y*x + j*x   + k] - V[5*3 + 1] * Ex[i*y*x + j*x   + k] + V[5*3 + 0] * Ey[i*y*x + j*x   + k]) - hz[i*y*x*q + j*x*q + k*q   + 5];
					hz[i*y*x*q + j*x*q + k*q  + 6] = 2 * W * (Hz[i*y*x + j*x   + k] - V[6*3 + 1] * Ex[i*y*x + j*x   + k] + V[6*3 + 0] * Ey[i*y*x + j*x   + k]) - hz[i*y*x*q + j*x*q + k*q   + 6];
				}
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


void collStream(float* ex, float* ey, float* ez, float* hx, float* hy, float* hz, float* Ex, float* Ey, float* Ez, float* Hx, float* Hy, float* Hz, float* er, float* mur, int z, int y, int x, int q) {

	int i=0, j=0, k=0;
	
	

	for (i = 0; i < z; i++) {
		for (j = 0; j < y; j++) {
			for (k = 0; k < x; k++) {
			
				int k1 = x-1-k, k11 = k1 - 1;
				int k21 = k+1;
				int j3 = y-1-j, j31 = j3-1;
				int j41 = j+1;

				if (k1 == 0)
					k11 = k1;
				if (k == x-1)
					k21 = k;
				if (j3 == 0)
					j31 = j3;
				if (j == y-1)
					j41 = j;
					
				ex[i*y*x*q + j*x*q + k*q  + 0] = 2 * ((er[i*y*x + j*x + k] - 1) * Ex[i*y*x + j*x + k]) - ex[i*y*x*q + j*x*q + k*q + 0];
				ex[i*y*x*q + j*x*q + k1*q + 1] = 2 * W * (Ex[i*y*x + j*x + k11] - V[1*3 + 1] * Hz[i*y*x + j*x + k11] + V[1*3 + 2] * Hy[i*y*x + j*x + k11]) - ex[i*y*x*q + j*x*q + k11*q + 1];
				ex[i*y*x*q + j*x*q + k*q  + 2] = 2 * W * (Ex[i*y*x + j*x + k21] - V[2*3 + 1] * Hz[i*y*x + j*x + k21] + V[2*3 + 2] * Hy[i*y*x + j*x + k21]) - ex[i*y*x*q + j*x*q + k21*q + 2];
				ex[i*y*x*q + j3*x*q + k*q + 3] = 2 * W * (Ex[i*y*x + j31*x + k] - V[3*3 + 1] * Hz[i*y*x + j31*x + k] + V[3*3 + 2] * Hy[i*y*x + j31*x + k]) - ex[i*y*x*q + j31*x*q + k*q + 3];
				ex[i*y*x*q + j*x*q + k*q  + 4] = 2 * W * (Ex[i*y*x + j41*x + k] - V[4*3 + 1] * Hz[i*y*x + j41*x + k] + V[4*3 + 2] * Hy[i*y*x + j41*x + k]) - ex[i*y*x*q + j41*x*q + k*q + 4];
				ex[i*y*x*q + j*x*q + k*q  + 5] = 2 * W * (Ex[i*y*x + j*x   + k] - V[5*3 + 1] * Hz[i*y*x + j*x   + k] + V[5*3 + 2] * Hy[i*y*x + j*x   + k]) - ex[i*y*x*q + j*x*q + k*q   + 5];
				ex[i*y*x*q + j*x*q + k*q  + 6] = 2 * W * (Ex[i*y*x + j*x   + k] - V[6*3 + 1] * Hz[i*y*x + j*x   + k] + V[6*3 + 2] * Hy[i*y*x + j*x   + k]) - ex[i*y*x*q + j*x*q + k*q   + 6];
				
				
				ey[i*y*x*q + j*x*q + k*q  + 0] = 2 * ((er[i*y*x + j*x + k] - 1) * Ey[i*y*x + j*x + k]) - ey[i*y*x*q + j*x*q + k*q + 0];
				ey[i*y*x*q + j*x*q + k1*q + 1] = 2 * W * (Ey[i*y*x + j*x + k11] - V[1*3 + 2] * Hx[i*y*x + j*x + k11] + V[1*3 + 0] * Hz[i*y*x + j*x + k11]) - ey[i*y*x*q + j*x*q + k11*q + 1];
				ey[i*y*x*q + j*x*q + k*q  + 2] = 2 * W * (Ey[i*y*x + j*x + k21] - V[2*3 + 2] * Hx[i*y*x + j*x + k21] + V[2*3 + 0] * Hz[i*y*x + j*x + k21]) - ey[i*y*x*q + j*x*q + k21*q + 2];
				ey[i*y*x*q + j3*x*q + k*q + 3] = 2 * W * (Ey[i*y*x + j31*x + k] - V[3*3 + 2] * Hx[i*y*x + j31*x + k] + V[3*3 + 0] * Hz[i*y*x + j31*x + k]) - ey[i*y*x*q + j31*x*q + k*q + 3];
				ey[i*y*x*q + j*x*q + k*q  + 4] = 2 * W * (Ey[i*y*x + j41*x + k] - V[4*3 + 2] * Hx[i*y*x + j41*x + k] + V[4*3 + 0] * Hz[i*y*x + j41*x + k]) - ey[i*y*x*q + j41*x*q + k*q + 4];
				ey[i*y*x*q + j*x*q + k*q  + 5] = 2 * W * (Ey[i*y*x + j*x   + k] - V[5*3 + 2] * Hx[i*y*x + j*x   + k] + V[5*3 + 0] * Hz[i*y*x + j*x   + k]) - ey[i*y*x*q + j*x*q + k*q   + 5];
				ey[i*y*x*q + j*x*q + k*q  + 6] = 2 * W * (Ey[i*y*x + j*x   + k] - V[6*3 + 2] * Hx[i*y*x + j*x   + k] + V[6*3 + 0] * Hz[i*y*x + j*x   + k]) - ey[i*y*x*q + j*x*q + k*q   + 6];


				ez[i*y*x*q + j*x*q + k*q  + 0] = 2 * ((er[i*y*x + j*x + k] - 1) * Ez[i*y*x + j*x + k]) - ez[i*y*x*q + j*x*q + k*q + 0];
				ez[i*y*x*q + j*x*q + k1*q + 1] = 2 * W * (Ez[i*y*x + j*x + k11] - V[1*3 + 0] * Hy[i*y*x + j*x + k11] + V[1*3 + 1] * Hx[i*y*x + j*x + k11]) - ez[i*y*x*q + j*x*q + k11*q + 1];
				ez[i*y*x*q + j*x*q + k*q  + 2] = 2 * W * (Ez[i*y*x + j*x + k21] - V[2*3 + 0] * Hy[i*y*x + j*x + k21] + V[2*3 + 1] * Hx[i*y*x + j*x + k21]) - ez[i*y*x*q + j*x*q + k21*q + 2];
				ez[i*y*x*q + j3*x*q + k*q + 3] = 2 * W * (Ez[i*y*x + j31*x + k] - V[3*3 + 0] * Hy[i*y*x + j31*x + k] + V[3*3 + 1] * Hx[i*y*x + j31*x + k]) - ez[i*y*x*q + j31*x*q + k*q + 3];
				ez[i*y*x*q + j*x*q + k*q  + 4] = 2 * W * (Ez[i*y*x + j41*x + k] - V[4*3 + 0] * Hy[i*y*x + j41*x + k] + V[4*3 + 1] * Hx[i*y*x + j41*x + k]) - ez[i*y*x*q + j41*x*q + k*q + 4];
				ez[i*y*x*q + j*x*q + k*q  + 5] = 2 * W * (Ez[i*y*x + j*x   + k] - V[5*3 + 0] * Hy[i*y*x + j*x   + k] + V[5*3 + 1] * Hx[i*y*x + j*x   + k]) - ez[i*y*x*q + j*x*q + k*q   + 5];
				ez[i*y*x*q + j*x*q + k*q  + 6] = 2 * W * (Ez[i*y*x + j*x   + k] - V[6*3 + 0] * Hy[i*y*x + j*x   + k] + V[6*3 + 1] * Hx[i*y*x + j*x   + k]) - ez[i*y*x*q + j*x*q + k*q   + 6];


				hx[i*y*x*q + j*x*q + k*q  + 0] = 2 * ((mur[i*y*x + j*x + k] - 1) * Hx[i*y*x + j*x + k]) - hx[i*y*x*q + j*x*q + k*q + 0];
				hx[i*y*x*q + j*x*q + k1*q + 1] = 2 * W * (Hx[i*y*x + j*x + k11] - V[1*3 + 2] * Ey[i*y*x + j*x + k11] + V[1*3 + 1] * Ez[i*y*x + j*x + k11]) - hx[i*y*x*q + j*x*q + k11*q + 1];
				hx[i*y*x*q + j*x*q + k*q  + 2] = 2 * W * (Hx[i*y*x + j*x + k21] - V[2*3 + 2] * Ey[i*y*x + j*x + k21] + V[2*3 + 1] * Ez[i*y*x + j*x + k21]) - hx[i*y*x*q + j*x*q + k21*q + 2];
				hx[i*y*x*q + j3*x*q + k*q + 3] = 2 * W * (Hx[i*y*x + j31*x + k] - V[3*3 + 2] * Ey[i*y*x + j31*x + k] + V[3*3 + 1] * Ez[i*y*x + j31*x + k]) - hx[i*y*x*q + j31*x*q + k*q + 3];
				hx[i*y*x*q + j*x*q + k*q  + 4] = 2 * W * (Hx[i*y*x + j41*x + k] - V[4*3 + 2] * Ey[i*y*x + j41*x + k] + V[4*3 + 1] * Ez[i*y*x + j41*x + k]) - hx[i*y*x*q + j41*x*q + k*q + 4];
				hx[i*y*x*q + j*x*q + k*q  + 5] = 2 * W * (Hx[i*y*x + j*x   + k] - V[5*3 + 2] * Ey[i*y*x + j*x   + k] + V[5*3 + 1] * Ez[i*y*x + j*x   + k]) - hx[i*y*x*q + j*x*q + k*q   + 5];
				hx[i*y*x*q + j*x*q + k*q  + 6] = 2 * W * (Hx[i*y*x + j*x   + k] - V[6*3 + 2] * Ey[i*y*x + j*x   + k] + V[6*3 + 1] * Ez[i*y*x + j*x   + k]) - hx[i*y*x*q + j*x*q + k*q   + 6];
				
				
				hy[i*y*x*q + j*x*q + k*q  + 0] = 2 * ((mur[i*y*x + j*x + k] - 1) * Hy[i*y*x + j*x + k]) - hy[i*y*x*q + j*x*q + k*q + 0];
				hy[i*y*x*q + j*x*q + k1*q + 1] = 2 * W * (Hy[i*y*x + j*x + k11] - V[1*3 + 0] * Ez[i*y*x + j*x + k11] + V[1*3 + 2] * Ex[i*y*x + j*x + k11]) - hy[i*y*x*q + j*x*q + k11*q + 1];
				hy[i*y*x*q + j*x*q + k*q  + 2] = 2 * W * (Hy[i*y*x + j*x + k21] - V[2*3 + 0] * Ez[i*y*x + j*x + k21] + V[2*3 + 2] * Ex[i*y*x + j*x + k21]) - hy[i*y*x*q + j*x*q + k21*q + 2];
				hy[i*y*x*q + j3*x*q + k*q + 3] = 2 * W * (Hy[i*y*x + j31*x + k] - V[3*3 + 0] * Ez[i*y*x + j31*x + k] + V[3*3 + 2] * Ex[i*y*x + j31*x + k]) - hy[i*y*x*q + j31*x*q + k*q + 3];
				hy[i*y*x*q + j*x*q + k*q  + 4] = 2 * W * (Hy[i*y*x + j41*x + k] - V[4*3 + 0] * Ez[i*y*x + j41*x + k] + V[4*3 + 2] * Ex[i*y*x + j41*x + k]) - hy[i*y*x*q + j41*x*q + k*q + 4];
				hy[i*y*x*q + j*x*q + k*q  + 5] = 2 * W * (Hy[i*y*x + j*x   + k] - V[5*3 + 0] * Ez[i*y*x + j*x   + k] + V[5*3 + 2] * Ex[i*y*x + j*x   + k]) - hy[i*y*x*q + j*x*q + k*q   + 5];
				hy[i*y*x*q + j*x*q + k*q  + 6] = 2 * W * (Hy[i*y*x + j*x   + k] - V[6*3 + 0] * Ez[i*y*x + j*x   + k] + V[6*3 + 2] * Ex[i*y*x + j*x   + k]) - hy[i*y*x*q + j*x*q + k*q   + 6];


				hz[i*y*x*q + j*x*q + k*q  + 0] = 2 * ((mur[i*y*x + j*x + k] - 1) * Hz[i*y*x + j*x + k]) - hz[i*y*x*q + j*x*q + k*q + 0];
				hz[i*y*x*q + j*x*q + k1*q + 1] = 2 * W * (Hz[i*y*x + j*x + k11] - V[1*3 + 1] * Ex[i*y*x + j*x + k11] + V[1*3 + 0] * Ey[i*y*x + j*x + k11]) - hz[i*y*x*q + j*x*q + k11*q + 1];
				hz[i*y*x*q + j*x*q + k*q  + 2] = 2 * W * (Hz[i*y*x + j*x + k21] - V[2*3 + 1] * Ex[i*y*x + j*x + k21] + V[2*3 + 0] * Ey[i*y*x + j*x + k21]) - hz[i*y*x*q + j*x*q + k21*q + 2];
				hz[i*y*x*q + j3*x*q + k*q + 3] = 2 * W * (Hz[i*y*x + j31*x + k] - V[3*3 + 1] * Ex[i*y*x + j31*x + k] + V[3*3 + 0] * Ey[i*y*x + j31*x + k]) - hz[i*y*x*q + j31*x*q + k*q + 3];
				hz[i*y*x*q + j*x*q + k*q  + 4] = 2 * W * (Hz[i*y*x + j41*x + k] - V[4*3 + 1] * Ex[i*y*x + j41*x + k] + V[4*3 + 0] * Ey[i*y*x + j41*x + k]) - hz[i*y*x*q + j41*x*q + k*q + 4];
				hz[i*y*x*q + j*x*q + k*q  + 5] = 2 * W * (Hz[i*y*x + j*x   + k] - V[5*3 + 1] * Ex[i*y*x + j*x   + k] + V[5*3 + 0] * Ey[i*y*x + j*x   + k]) - hz[i*y*x*q + j*x*q + k*q   + 5];
				hz[i*y*x*q + j*x*q + k*q  + 6] = 2 * W * (Hz[i*y*x + j*x   + k] - V[6*3 + 1] * Ex[i*y*x + j*x   + k] + V[6*3 + 0] * Ey[i*y*x + j*x   + k]) - hz[i*y*x*q + j*x*q + k*q   + 6];
			}
		}
	}
}


/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/





