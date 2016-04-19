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
	
	/*Initiation of grid boundary */
	for(i = 0; i < gridSize; i++){
		grid[0][i] = 1;				// upper boundary
		new[0][i] = 1;
		grid[gridSize-1][i] = 1;	// lower boundary
		new[gridSize-1][i] = 1;
		grid[i][0] = 1; 			// left  boundary
		new[i][0] = 1;
		grid[i][gridSize-1] = 1;	// right boundary
		new[i][gridSize-1] = 1;
	}
	
	/*Instantiation of grids */
	for(i = 1; i < gridSize-1; i++){
		for(j = 1; j < gridSize-1; j++){
			grid[i][j] = 0;
		}	
	}
	
	/*Rows that will be processed by each processor*/
	
	workLoad = gridSize / numProc;
	
	/*Vi bör eventuellt använda shift istället för * 0.25*/
	/*Jacobi iterations*/
	
	for(k = 0; k < numIters; k++) {
		
		#pragma omp parallel for schedule(static, workLoad) private(j)
		for(i = 1; i < gridSize-1; i++) {
			for(j = 1; j < gridSize-1; j++) {
				new[i][j] = (grid[i-1][j] +
							 grid[i+1][j] +
							 grid[i][j-1] +
							 grid[i][j+1]) * 0.25;
			}
		}
		/*implicit barrier (it is needed)*/
		#pragma omp parallel for schedule(static, workLoad) private(j)
		for(i = 1; i < gridSize-1; i++) {
			for(j = 1; j < gridSize-1; j++) {
				grid[i][j] = (new[i-1][j] +
							  new[i+1][j] +
							  new[i][j-1] +
							  new[i][j+1]) * 0.25;
			}
		}
		/*implicit barrier*/
	}
	
	/*computation of the maximum difference*/
	maxdiff = 0;
	for(i = 1; i < gridSize-1; i++){
		for(j = 1; j < gridSize-1; j++){
			temp = grid[i][j] - new[i][j];
			if(temp < 0)
				temp = -temp;
			if(temp > maxdiff)
				maxdiff = temp;
		}
	}
	/*print the maximum difference*/
	printf("Maximum difference: %f\n", maxdiff);
	
	/*Prints the grid with boundary*/
	for(i = 0; i < gridSize; i++){
			printf("\n");
		for(j = 0; j < gridSize; j++){
			printf("%f  ", grid[i][j]);
		}
	}
	
	/*Behöver vi kanske köra free på vartenda element eller duger detta tror du?*/
	free(grid);
	free(new);
	
	printf("\n");
	exit(0);
}
