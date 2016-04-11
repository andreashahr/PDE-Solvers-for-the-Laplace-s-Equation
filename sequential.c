/*
	Author Kardå Chalak & Andreas Hahr
	kardac@kth.se
	andreas epost
	
	A sequential Jacobi iteration program
*/
 #include <stdlib.h>
 #include <stdio.h>
 
 #define GRIDSIZE 10
 #define NUMITERS 10
 
 int main(int argc, char *argv[]){
 
	/*gridSize, the grid size not including boundaries
	  numIters, the number of iterations to use
	 */
	
	int gridSize, numIters, i, j, k;
	double temp, maxdiff;
	
	gridSize = (argc > 1)? atoi(argv[1]) + 2 : GRIDSIZE;
	numIters = (argc > 2)? atoi(argv[2]) : NUMITERS;
	
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
	
	/*Vi bör eventuellt använda shift istället för * 0.25*/
	/*Jacobi iterations*/
	for(k = 0; k < numIters; k++) {
	
		for(i = 1; i < gridSize-1; i++) {
			for(j = 1; j < gridSize-1; j++) {
				new[i][j] = (grid[i-1][j] +
							 grid[i+1][j] +
							 grid[i][j-1] +
							 grid[i][j+1]) * 0.25;
			}
		}
		for(i = 1; i < gridSize-1; i++) {
			for(j = 1; j < gridSize-1; j++) {
				grid[i][j] = (new[i-1][j] +
							  new[i+1][j] +
							  new[i][j-1] +
							  new[i][j+1]) * 0.25;
			}
		}
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
