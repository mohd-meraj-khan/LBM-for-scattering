#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define W 1.0/6.0



int V[21] = {0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1};


/*---------------------------------------------------------------------------------------------*/
// calculates the macroscopic fields from the distribution functions


void macroField(float* dis_func, float* mat_prop, float* field, int y, int x, int q) {
	
	int i, j, k;
	
	for(i = 0; i < y; i++) {
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


void collStreamForcingNode(float* ex, float* ey, float* ez, float* hx, float* hy, float* hz, float* Ex, float* Ey, float* Ez, float* Hx, float* Hy, float* Hz, float* er, float* mur, int y, int x, int q, int xloc, int ymin, int ymax) {

	int i, j;
	
	
	for (i = 0; i < y; i++) {
		for (j = 0; j < x; j++) {
		
			int j1 = x-1-j, j11 = j1-1;
			int j21 = j+1;
			int i3 = y-1-i, i31 = i3-1;
			int i41 = i+1;

			if (j1 == 0)
				j11 = j1;
			if (j == x-1)
				j21 = j;
			if (i3 == 0)
				i31 = i3;
			if (i == y-1)
				i41 = i;
			
			if (i >= ymin && i < ymax && j == xloc) { // no collision at the forcing node

				ex[i*x*q + j*q  + 0] = ((er[i*x + j] - 1) * Ex[i*x + j]);	
				ex[i*x*q + j1*q + 1] = W * (Ex[i*x + j11] - V[1*3 + 1] * Hz[i*x + j11] + V[1*3 + 2] * Hy[i*x + j11]);
				ex[i*x*q + j*q  + 2] = W * (Ex[i*x + j21] - V[2*3 + 1] * Hz[i*x + j21] + V[2*3 + 2] * Hy[i*x + j21]);
				ex[i3*x*q + j*q + 3] = W * (Ex[i31*x + j] - V[3*3 + 1] * Hz[i31*x + j] + V[3*3 + 2] * Hy[i31*x + j]);
				ex[i*x*q + j*q  + 4] = W * (Ex[i41*x + j] - V[4*3 + 1] * Hz[i41*x + j] + V[4*3 + 2] * Hy[i41*x + j]);
				ex[i*x*q + j*q  + 5] = W * (Ex[i*x   + j] - V[5*3 + 1] * Hz[i*x   + j] + V[5*3 + 2] * Hy[i*x   + j]);
				ex[i*x*q + j*q  + 6] = W * (Ex[i*x   + j] - V[6*3 + 1] * Hz[i*x   + j] + V[6*3 + 2] * Hy[i*x   + j]);
				
				
				ey[i*x*q + j*q  + 0] = ((er[i*x + j] - 1) * Ey[i*x + j]);
				ey[i*x*q + j1*q + 1] = W * (Ey[i*x + j11] - V[1*3 + 2] * Hx[i*x + j11] + V[1*3 + 0] * Hz[i*x + j11]);
				ey[i*x*q + j*q  + 2] = W * (Ey[i*x + j21] - V[2*3 + 2] * Hx[i*x + j21] + V[2*3 + 0] * Hz[i*x + j21]);
				ey[i3*x*q + j*q + 3] = W * (Ey[i31*x + j] - V[3*3 + 2] * Hx[i31*x + j] + V[3*3 + 0] * Hz[i31*x + j]);
				ey[i*x*q + j*q  + 4] = W * (Ey[i41*x + j] - V[4*3 + 2] * Hx[i41*x + j] + V[4*3 + 0] * Hz[i41*x + j]);
				ey[i*x*q + j*q  + 5] = W * (Ey[i*x   + j] - V[5*3 + 2] * Hx[i*x   + j] + V[5*3 + 0] * Hz[i*x   + j]);
				ey[i*x*q + j*q  + 6] = W * (Ey[i*x   + j] - V[6*3 + 2] * Hx[i*x   + j] + V[6*3 + 0] * Hz[i*x   + j]);


				ez[i*x*q + j*q  + 0] = ((er[i*x + j] - 1) * Ez[i*x + j]);
				ez[i*x*q + j1*q + 1] = W * (Ez[i*x + j11] - V[1*3 + 0] * Hy[i*x + j11] + V[1*3 + 1] * Hx[i*x + j11]);
				ez[i*x*q + j*q  + 2] = W * (Ez[i*x + j21] - V[2*3 + 0] * Hy[i*x + j21] + V[2*3 + 1] * Hx[i*x + j21]);
				ez[i3*x*q + j*q + 3] = W * (Ez[i31*x + j] - V[3*3 + 0] * Hy[i31*x + j] + V[3*3 + 1] * Hx[i31*x + j]);
				ez[i*x*q + j*q  + 4] = W * (Ez[i41*x + j] - V[4*3 + 0] * Hy[i41*x + j] + V[4*3 + 1] * Hx[i41*x + j]);
				ez[i*x*q + j*q  + 5] = W * (Ez[i*x + j]   - V[5*3 + 0] * Hy[i*x   + j] + V[5*3 + 1] * Hx[i*x   + j]);
				ez[i*x*q + j*q  + 6] = W * (Ez[i*x + j]   - V[6*3 + 0] * Hy[i*x   + j] + V[6*3 + 1] * Hx[i*x   + j]);



				hx[i*x*q + j*q  + 0] = ((mur[i*x + j] - 1) * Hx[i*x + j]);
				hx[i*x*q + j1*q + 1] = W * (Hx[i*x + j11] - V[1*3 + 2] * Ey[i*x + j11] + V[1*3 + 1] * Ez[i*x + j11]);
				hx[i*x*q + j*q  + 2] = W * (Hx[i*x + j21] - V[2*3 + 2] * Ey[i*x + j21] + V[2*3 + 1] * Ez[i*x + j21]);
				hx[i3*x*q + j*q + 3] = W * (Hx[i31*x + j] - V[3*3 + 2] * Ey[i31*x + j] + V[3*3 + 1] * Ez[i31*x + j]);
				hx[i*x*q + j*q  + 4] = W * (Hx[i41*x + j] - V[4*3 + 2] * Ey[i41*x + j] + V[4*3 + 1] * Ez[i41*x + j]);
				hx[i*x*q + j*q  + 5] = W * (Hx[i*x   + j] - V[5*3 + 2] * Ey[i*x   + j] + V[5*3 + 1] * Ez[i*x   + j]);
				hx[i*x*q + j*q  + 6] = W * (Hx[i*x   + j] - V[6*3 + 2] * Ey[i*x   + j] + V[6*3 + 1] * Ez[i*x   + j]);
				
				
			
				hy[i*x*q + j*q  + 0] = ((mur[i*x + j] - 1) * Hy[i*x + j]);
				hy[i*x*q + j1*q + 1] = W * (Hy[i*x + j11] - V[1*3 + 0] * Ez[i*x + j11] + V[1*3 + 2] * Ex[i*x + j11]);
				hy[i*x*q + j*q  + 2] = W * (Hy[i*x + j21] - V[2*3 + 0] * Ez[i*x + j21] + V[2*3 + 2] * Ex[i*x + j21]);
				hy[i3*x*q + j*q + 3] = W * (Hy[i31*x + j] - V[3*3 + 0] * Ez[i31*x + j] + V[3*3 + 2] * Ex[i31*x + j]);
				hy[i*x*q + j*q  + 4] = W * (Hy[i41*x + j] - V[4*3 + 0] * Ez[i41*x + j] + V[4*3 + 2] * Ex[i41*x + j]);
				hy[i*x*q + j*q  + 5] = W * (Hy[i*x   + j] - V[5*3 + 0] * Ez[i*x   + j] + V[5*3 + 2] * Ex[i*x   + j]);
				hy[i*x*q + j*q  + 6] = W * (Hy[i*x   + j] - V[6*3 + 0] * Ez[i*x   + j] + V[6*3 + 2] * Ex[i*x   + j]);


				hz[i*x*q + j*q  + 0] = ((mur[i*x + j] - 1) * Hz[i*x + j]);
				hz[i*x*q + j1*q + 1] = W * (Hz[i*x + j11] - V[1*3 + 1] * Ex[i*x + j11] + V[1*3 + 0] * Ey[i*x + j11]);
				hz[i*x*q + j*q  + 2] = W * (Hz[i*x + j21] - V[2*3 + 1] * Ex[i*x + j21] + V[2*3 + 0] * Ey[i*x + j21]);
				hz[i3*x*q + j*q + 3] = W * (Hz[i31*x + j] - V[3*3 + 1] * Ex[i31*x + j] + V[3*3 + 0] * Ey[i31*x + j]);
				hz[i*x*q + j*q  + 4] = W * (Hz[i41*x + j] - V[4*3 + 1] * Ex[i41*x + j] + V[4*3 + 0] * Ey[i41*x + j]);
				hz[i*x*q + j*q  + 5] = W * (Hz[i*x   + j] - V[5*3 + 1] * Ex[i*x   + j] + V[5*3 + 0] * Ey[i*x   + j]);
				hz[i*x*q + j*q  + 6] = W * (Hz[i*x   + j] - V[6*3 + 1] * Ex[i*x   + j] + V[6*3 + 0] * Ey[i*x   + j]);
			}
			
							
			
			else {
				
				ex[i*x*q + j*q  + 0] = 2 * ((er[i*x + j] - 1) * Ex[i*x + j]) - ex[i*x*q + j*q + 0];
				ex[i*x*q + j1*q + 1] = 2 * W * (Ex[i*x + j11] - V[1*3 + 1] * Hz[i*x + j11] + V[1*3 + 2] * Hy[i*x + j11]) - ex[i*x*q + j11*q + 1];
				ex[i*x*q + j*q  + 2] = 2 * W * (Ex[i*x + j21] - V[2*3 + 1] * Hz[i*x + j21] + V[2*3 + 2] * Hy[i*x + j21]) - ex[i*x*q + j21*q + 2];
				ex[i3*x*q + j*q + 3] = 2 * W * (Ex[i31*x + j] - V[3*3 + 1] * Hz[i31*x + j] + V[3*3 + 2] * Hy[i31*x + j]) - ex[i31*x*q + j*q + 3];
				ex[i*x*q + j*q  + 4] = 2 * W * (Ex[i41*x + j] - V[4*3 + 1] * Hz[i41*x + j] + V[4*3 + 2] * Hy[i41*x + j]) - ex[i41*x*q + j*q + 4];
				ex[i*x*q + j*q  + 5] = 2 * W * (Ex[i*x   + j] - V[5*3 + 1] * Hz[i*x   + j] + V[5*3 + 2] * Hy[i*x   + j]) - ex[i*x*q + j*q   + 5];
				ex[i*x*q + j*q  + 6] = 2 * W * (Ex[i*x   + j] - V[6*3 + 1] * Hz[i*x   + j] + V[6*3 + 2] * Hy[i*x   + j]) - ex[i*x*q + j*q   + 6];
				
				
				ey[i*x*q + j*q  + 0] = 2 * ((er[i*x + j] - 1) * Ey[i*x + j]) - ey[i*x*q + j*q + 0];
				ey[i*x*q + j1*q + 1] = 2 * W * (Ey[i*x + j11] - V[1*3 + 2] * Hx[i*x + j11] + V[1*3 + 0] * Hz[i*x + j11]) - ey[i*x*q + j11*q + 1];
				ey[i*x*q + j*q  + 2] = 2 * W * (Ey[i*x + j21] - V[2*3 + 2] * Hx[i*x + j21] + V[2*3 + 0] * Hz[i*x + j21]) - ey[i*x*q + j21*q + 2];
				ey[i3*x*q + j*q + 3] = 2 * W * (Ey[i31*x + j] - V[3*3 + 2] * Hx[i31*x + j] + V[3*3 + 0] * Hz[i31*x + j]) - ey[i31*x*q + j*q + 3];
				ey[i*x*q + j*q  + 4] = 2 * W * (Ey[i41*x + j] - V[4*3 + 2] * Hx[i41*x + j] + V[4*3 + 0] * Hz[i41*x + j]) - ey[i41*x*q + j*q + 4];
				ey[i*x*q + j*q  + 5] = 2 * W * (Ey[i*x   + j] - V[5*3 + 2] * Hx[i*x   + j] + V[5*3 + 0] * Hz[i*x   + j]) - ey[i*x*q + j*q   + 5];
				ey[i*x*q + j*q  + 6] = 2 * W * (Ey[i*x   + j] - V[6*3 + 2] * Hx[i*x   + j] + V[6*3 + 0] * Hz[i*x   + j]) - ey[i*x*q + j*q   + 6];


				ez[i*x*q + j*q  + 0] = 2 * ((er[i*x + j] - 1) * Ez[i*x + j]) - ez[i*x*q + j*q + 0];
				ez[i*x*q + j1*q + 1] = 2 * W * (Ez[i*x + j11] - V[1*3 + 0] * Hy[i*x + j11] + V[1*3 + 1] * Hx[i*x + j11]) - ez[i*x*q + j11*q + 1];
				ez[i*x*q + j*q  + 2] = 2 * W * (Ez[i*x + j21] - V[2*3 + 0] * Hy[i*x + j21] + V[2*3 + 1] * Hx[i*x + j21]) - ez[i*x*q + j21*q + 2];
				ez[i3*x*q + j*q + 3] = 2 * W * (Ez[i31*x + j] - V[3*3 + 0] * Hy[i31*x + j] + V[3*3 + 1] * Hx[i31*x + j]) - ez[i31*x*q + j*q + 3];
				ez[i*x*q + j*q  + 4] = 2 * W * (Ez[i41*x + j] - V[4*3 + 0] * Hy[i41*x + j] + V[4*3 + 1] * Hx[i41*x + j]) - ez[i41*x*q + j*q + 4];
				ez[i*x*q + j*q  + 5] = 2 * W * (Ez[i*x   + j] - V[5*3 + 0] * Hy[i*x   + j] + V[5*3 + 1] * Hx[i*x   + j]) - ez[i*x*q + j*q   + 5];
				ez[i*x*q + j*q  + 6] = 2 * W * (Ez[i*x   + j] - V[6*3 + 0] * Hy[i*x   + j] + V[6*3 + 1] * Hx[i*x   + j]) - ez[i*x*q + j*q   + 6];


				hx[i*x*q + j*q  + 0] = 2 * ((mur[i*x + j] - 1) * Hx[i*x + j]) - hx[i*x*q + j*q + 0];
				hx[i*x*q + j1*q + 1] = 2 * W * (Hx[i*x + j11] - V[1*3 + 2] * Ey[i*x + j11] + V[1*3 + 1] * Ez[i*x + j11]) - hx[i*x*q + j11*q + 1];
				hx[i*x*q + j*q  + 2] = 2 * W * (Hx[i*x + j21] - V[2*3 + 2] * Ey[i*x + j21] + V[2*3 + 1] * Ez[i*x + j21]) - hx[i*x*q + j21*q + 2];
				hx[i3*x*q + j*q + 3] = 2 * W * (Hx[i31*x + j] - V[3*3 + 2] * Ey[i31*x + j] + V[3*3 + 1] * Ez[i31*x + j]) - hx[i31*x*q + j*q + 3];
				hx[i*x*q + j*q  + 4] = 2 * W * (Hx[i41*x + j] - V[4*3 + 2] * Ey[i41*x + j] + V[4*3 + 1] * Ez[i41*x + j]) - hx[i41*x*q + j*q + 4];
				hx[i*x*q + j*q  + 5] = 2 * W * (Hx[i*x   + j] - V[5*3 + 2] * Ey[i*x   + j] + V[5*3 + 1] * Ez[i*x   + j]) - hx[i*x*q + j*q   + 5];
				hx[i*x*q + j*q  + 6] = 2 * W * (Hx[i*x   + j] - V[6*3 + 2] * Ey[i*x   + j] + V[6*3 + 1] * Ez[i*x   + j]) - hx[i*x*q + j*q   + 6];
				
				
				hy[i*x*q + j*q  + 0] = 2 * ((mur[i*x + j] - 1) * Hy[i*x + j]) - hy[i*x*q + j*q + 0];
				hy[i*x*q + j1*q + 1] = 2 * W * (Hy[i*x + j11] - V[1*3 + 0] * Ez[i*x + j11] + V[1*3 + 2] * Ex[i*x + j11]) - hy[i*x*q + j11*q + 1];
				hy[i*x*q + j*q  + 2] = 2 * W * (Hy[i*x + j21] - V[2*3 + 0] * Ez[i*x + j21] + V[2*3 + 2] * Ex[i*x + j21]) - hy[i*x*q + j21*q + 2];
				hy[i3*x*q + j*q + 3] = 2 * W * (Hy[i31*x + j] - V[3*3 + 0] * Ez[i31*x + j] + V[3*3 + 2] * Ex[i31*x + j]) - hy[i31*x*q + j*q + 3];
				hy[i*x*q + j*q  + 4] = 2 * W * (Hy[i41*x + j] - V[4*3 + 0] * Ez[i41*x + j] + V[4*3 + 2] * Ex[i41*x + j]) - hy[i41*x*q + j*q + 4];
				hy[i*x*q + j*q  + 5] = 2 * W * (Hy[i*x   + j] - V[5*3 + 0] * Ez[i*x   + j] + V[5*3 + 2] * Ex[i*x   + j]) - hy[i*x*q + j*q   + 5];
				hy[i*x*q + j*q  + 6] = 2 * W * (Hy[i*x   + j] - V[6*3 + 0] * Ez[i*x   + j] + V[6*3 + 2] * Ex[i*x   + j]) - hy[i*x*q + j*q   + 6];


				hz[i*x*q + j*q  + 0] = 2 * ((mur[i*x + j] - 1) * Hz[i*x + j]) - hz[i*x*q + j*q + 0];
				hz[i*x*q + j1*q + 1] = 2 * W * (Hz[i*x + j11] - V[1*3 + 1] * Ex[i*x + j11] + V[1*3 + 0] * Ey[i*x + j11]) - hz[i*x*q + j11*q + 1];
				hz[i*x*q + j*q  + 2] = 2 * W * (Hz[i*x + j21] - V[2*3 + 1] * Ex[i*x + j21] + V[2*3 + 0] * Ey[i*x + j21]) - hz[i*x*q + j21*q + 2];
				hz[i3*x*q + j*q + 3] = 2 * W * (Hz[i31*x + j] - V[3*3 + 1] * Ex[i31*x + j] + V[3*3 + 0] * Ey[i31*x + j]) - hz[i31*x*q + j*q + 3];
				hz[i*x*q + j*q  + 4] = 2 * W * (Hz[i41*x + j] - V[4*3 + 1] * Ex[i41*x + j] + V[4*3 + 0] * Ey[i41*x + j]) - hz[i41*x*q + j*q + 4];
				hz[i*x*q + j*q  + 5] = 2 * W * (Hz[i*x   + j] - V[5*3 + 1] * Ex[i*x   + j] + V[5*3 + 0] * Ey[i*x   + j]) - hz[i*x*q + j*q   + 5];
				hz[i*x*q + j*q  + 6] = 2 * W * (Hz[i*x   + j] - V[6*3 + 1] * Ex[i*x   + j] + V[6*3 + 0] * Ey[i*x   + j]) - hz[i*x*q + j*q   + 6];
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


void collStream(float* ex, float* ey, float* ez, float* hx, float* hy, float* hz, float* Ex, float* Ey, float* Ez, float* Hx, float* Hy, float* Hz, float* er, float* mur, int y, int x, int q) {

	int i, j;
	
	

	
	for (i = 0; i < y; i++) {
		for (j = 0; j < x; j++) {
		
			int j1 = x-1-j, j11 = j1 - 1;
			int j21 = j+1;
			int i3 = y-1-i, i31 = i3-1;
			int i41 = i+1;

			if (j1 == 0)
				j11 = j1;
			if (j == x-1)
				j21 = j;
			if (i3 == 0)
				i31 = i3;
			if (i == y-1)
				i41 = i;
				
			ex[i*x*q + j*q  + 0] = 2 * ((er[i*x + j] - 1) * Ex[i*x + j]) - ex[i*x*q + j*q + 0];
			ex[i*x*q + j1*q + 1] = 2 * W * (Ex[i*x + j11] - V[1*3 + 1] * Hz[i*x + j11] + V[1*3 + 2] * Hy[i*x + j11]) - ex[i*x*q + j11*q + 1];
			ex[i*x*q + j*q  + 2] = 2 * W * (Ex[i*x + j21] - V[2*3 + 1] * Hz[i*x + j21] + V[2*3 + 2] * Hy[i*x + j21]) - ex[i*x*q + j21*q + 2];
			ex[i3*x*q + j*q + 3] = 2 * W * (Ex[i31*x + j] - V[3*3 + 1] * Hz[i31*x + j] + V[3*3 + 2] * Hy[i31*x + j]) - ex[i31*x*q + j*q + 3];
			ex[i*x*q + j*q  + 4] = 2 * W * (Ex[i41*x + j] - V[4*3 + 1] * Hz[i41*x + j] + V[4*3 + 2] * Hy[i41*x + j]) - ex[i41*x*q + j*q + 4];
			ex[i*x*q + j*q  + 5] = 2 * W * (Ex[i*x   + j] - V[5*3 + 1] * Hz[i*x   + j] + V[5*3 + 2] * Hy[i*x   + j]) - ex[i*x*q + j*q   + 5];
			ex[i*x*q + j*q  + 6] = 2 * W * (Ex[i*x   + j] - V[6*3 + 1] * Hz[i*x   + j] + V[6*3 + 2] * Hy[i*x   + j]) - ex[i*x*q + j*q   + 6];
			
			
			ey[i*x*q + j*q  + 0] = 2 * ((er[i*x + j] - 1) * Ey[i*x + j]) - ey[i*x*q + j*q + 0];
			ey[i*x*q + j1*q + 1] = 2 * W * (Ey[i*x + j11] - V[1*3 + 2] * Hx[i*x + j11] + V[1*3 + 0] * Hz[i*x + j11]) - ey[i*x*q + j11*q + 1];
			ey[i*x*q + j*q  + 2] = 2 * W * (Ey[i*x + j21] - V[2*3 + 2] * Hx[i*x + j21] + V[2*3 + 0] * Hz[i*x + j21]) - ey[i*x*q + j21*q + 2];
			ey[i3*x*q + j*q + 3] = 2 * W * (Ey[i31*x + j] - V[3*3 + 2] * Hx[i31*x + j] + V[3*3 + 0] * Hz[i31*x + j]) - ey[i31*x*q + j*q + 3];
			ey[i*x*q + j*q  + 4] = 2 * W * (Ey[i41*x + j] - V[4*3 + 2] * Hx[i41*x + j] + V[4*3 + 0] * Hz[i41*x + j]) - ey[i41*x*q + j*q + 4];
			ey[i*x*q + j*q  + 5] = 2 * W * (Ey[i*x   + j] - V[5*3 + 2] * Hx[i*x   + j] + V[5*3 + 0] * Hz[i*x   + j]) - ey[i*x*q + j*q   + 5];
			ey[i*x*q + j*q  + 6] = 2 * W * (Ey[i*x   + j] - V[6*3 + 2] * Hx[i*x   + j] + V[6*3 + 0] * Hz[i*x   + j]) - ey[i*x*q + j*q   + 6];


			ez[i*x*q + j*q  + 0] = 2 * ((er[i*x + j] - 1) * Ez[i*x + j]) - ez[i*x*q + j*q + 0];
			ez[i*x*q + j1*q + 1] = 2 * W * (Ez[i*x + j11] - V[1*3 + 0] * Hy[i*x + j11] + V[1*3 + 1] * Hx[i*x + j11]) - ez[i*x*q + j11*q + 1];
			ez[i*x*q + j*q  + 2] = 2 * W * (Ez[i*x + j21] - V[2*3 + 0] * Hy[i*x + j21] + V[2*3 + 1] * Hx[i*x + j21]) - ez[i*x*q + j21*q + 2];
			ez[i3*x*q + j*q + 3] = 2 * W * (Ez[i31*x + j] - V[3*3 + 0] * Hy[i31*x + j] + V[3*3 + 1] * Hx[i31*x + j]) - ez[i31*x*q + j*q + 3];
			ez[i*x*q + j*q  + 4] = 2 * W * (Ez[i41*x + j] - V[4*3 + 0] * Hy[i41*x + j] + V[4*3 + 1] * Hx[i41*x + j]) - ez[i41*x*q + j*q + 4];
			ez[i*x*q + j*q  + 5] = 2 * W * (Ez[i*x   + j] - V[5*3 + 0] * Hy[i*x   + j] + V[5*3 + 1] * Hx[i*x   + j]) - ez[i*x*q + j*q   + 5];
			ez[i*x*q + j*q  + 6] = 2 * W * (Ez[i*x   + j] - V[6*3 + 0] * Hy[i*x   + j] + V[6*3 + 1] * Hx[i*x   + j]) - ez[i*x*q + j*q   + 6];


			hx[i*x*q + j*q  + 0] = 2 * ((mur[i*x + j] - 1) * Hx[i*x + j]) - hx[i*x*q + j*q + 0];
			hx[i*x*q + j1*q + 1] = 2 * W * (Hx[i*x + j11] - V[1*3 + 2] * Ey[i*x + j11] + V[1*3 + 1] * Ez[i*x + j11]) - hx[i*x*q + j11*q + 1];
			hx[i*x*q + j*q  + 2] = 2 * W * (Hx[i*x + j21] - V[2*3 + 2] * Ey[i*x + j21] + V[2*3 + 1] * Ez[i*x + j21]) - hx[i*x*q + j21*q + 2];
			hx[i3*x*q + j*q + 3] = 2 * W * (Hx[i31*x + j] - V[3*3 + 2] * Ey[i31*x + j] + V[3*3 + 1] * Ez[i31*x + j]) - hx[i31*x*q + j*q + 3];
			hx[i*x*q + j*q  + 4] = 2 * W * (Hx[i41*x + j] - V[4*3 + 2] * Ey[i41*x + j] + V[4*3 + 1] * Ez[i41*x + j]) - hx[i41*x*q + j*q + 4];
			hx[i*x*q + j*q  + 5] = 2 * W * (Hx[i*x   + j] - V[5*3 + 2] * Ey[i*x   + j] + V[5*3 + 1] * Ez[i*x   + j]) - hx[i*x*q + j*q   + 5];
			hx[i*x*q + j*q  + 6] = 2 * W * (Hx[i*x   + j] - V[6*3 + 2] * Ey[i*x   + j] + V[6*3 + 1] * Ez[i*x   + j]) - hx[i*x*q + j*q   + 6];
			
			
			hy[i*x*q + j*q  + 0] = 2 * ((mur[i*x + j] - 1) * Hy[i*x + j]) - hy[i*x*q + j*q + 0];
			hy[i*x*q + j1*q + 1] = 2 * W * (Hy[i*x + j11] - V[1*3 + 0] * Ez[i*x + j11] + V[1*3 + 2] * Ex[i*x + j11]) - hy[i*x*q + j11*q + 1];
			hy[i*x*q + j*q  + 2] = 2 * W * (Hy[i*x + j21] - V[2*3 + 0] * Ez[i*x + j21] + V[2*3 + 2] * Ex[i*x + j21]) - hy[i*x*q + j21*q + 2];
			hy[i3*x*q + j*q + 3] = 2 * W * (Hy[i31*x + j] - V[3*3 + 0] * Ez[i31*x + j] + V[3*3 + 2] * Ex[i31*x + j]) - hy[i31*x*q + j*q + 3];
			hy[i*x*q + j*q  + 4] = 2 * W * (Hy[i41*x + j] - V[4*3 + 0] * Ez[i41*x + j] + V[4*3 + 2] * Ex[i41*x + j]) - hy[i41*x*q + j*q + 4];
			hy[i*x*q + j*q  + 5] = 2 * W * (Hy[i*x   + j] - V[5*3 + 0] * Ez[i*x   + j] + V[5*3 + 2] * Ex[i*x   + j]) - hy[i*x*q + j*q   + 5];
			hy[i*x*q + j*q  + 6] = 2 * W * (Hy[i*x   + j] - V[6*3 + 0] * Ez[i*x   + j] + V[6*3 + 2] * Ex[i*x   + j]) - hy[i*x*q + j*q   + 6];


			hz[i*x*q + j*q  + 0] = 2 * ((mur[i*x + j] - 1) * Hz[i*x + j]) - hz[i*x*q + j*q + 0];
			hz[i*x*q + j1*q + 1] = 2 * W * (Hz[i*x + j11] - V[1*3 + 1] * Ex[i*x + j11] + V[1*3 + 0] * Ey[i*x + j11]) - hz[i*x*q + j11*q + 1];
			hz[i*x*q + j*q  + 2] = 2 * W * (Hz[i*x + j21] - V[2*3 + 1] * Ex[i*x + j21] + V[2*3 + 0] * Ey[i*x + j21]) - hz[i*x*q + j21*q + 2];
			hz[i3*x*q + j*q + 3] = 2 * W * (Hz[i31*x + j] - V[3*3 + 1] * Ex[i31*x + j] + V[3*3 + 0] * Ey[i31*x + j]) - hz[i31*x*q + j*q + 3];
			hz[i*x*q + j*q  + 4] = 2 * W * (Hz[i41*x + j] - V[4*3 + 1] * Ex[i41*x + j] + V[4*3 + 0] * Ey[i41*x + j]) - hz[i41*x*q + j*q + 4];
			hz[i*x*q + j*q  + 5] = 2 * W * (Hz[i*x   + j] - V[5*3 + 1] * Ex[i*x   + j] + V[5*3 + 0] * Ey[i*x   + j]) - hz[i*x*q + j*q   + 5];
			hz[i*x*q + j*q  + 6] = 2 * W * (Hz[i*x   + j] - V[6*3 + 1] * Ex[i*x   + j] + V[6*3 + 0] * Ey[i*x   + j]) - hz[i*x*q + j*q   + 6];
		}
	}
}


/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/





