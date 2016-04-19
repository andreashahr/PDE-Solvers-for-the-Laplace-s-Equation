/*
	Author Kardå Chalak & Andreas Hahr
	kardac@kth.se
	andreas epost
	
	A parallel Jacobi iteration program using shared variables.
	Adapted for hardware with shared memory.
	
	Input: gridSize, numIters, numProc
	
	At the moment the program assumes that the number of rows is evenly
	divided by the number of processors. Gives slower performance when not true.
	Should be fixed!
*/

 #include <stdlib.h>
 #include <stdio.h>
 
 #include <omp.h>
 
 #define GRIDSIZE 10
 #define NUMITERS 10

int main(int argc, char *argv[]){

	/*gridSize, the grid size not including boundaries
	  numIters, the number of iterations to use
	*/
	
	int gridSize, numIters, numProc, workLoad, i, j, k;
	double temp, maxdiff;
	
	gridSize = (argc > 1)? atoi(argv[1]) + 2 : GRIDSIZE;
	numIters = (argc > 2)? atoi(argv[2]) : NUMITERS;
	
	/* sets threads*/
	int maxthreads = omp_get_max_threads();
	numProc = (argc > 2)? atoi(argv[2]) : maxthreads;
	if (numProc > maxthreads) numProc = maxthreads;
	if (numProc > gridSize) numProc = gridSize;
	
	printf("gridSize: %d\n", gridSize);
	printf("numIters: %d\n", numIters);
	
	/*malloc the columns*/
	double **grid = (double **) malloc(gridSize * sizeof(double*));
	double **new = (double **) malloc(gridSize * sizeof(double*));
	/*malloc the rows*/
	for(i = 0; i < gridSize; i++){
		grid[i] = (double *) malloc(gridSize * sizeof(double));
		new[i] = (double *) malloc(gridSize * sizeof(double));
	}
	
	void barrier(int id){
		// stub function for barrier
	}
	
	
	
	
	/*Behöver vi kanske köra free på vartenda element eller duger detta tror du?*/
	free(grid);
	free(new);
}
