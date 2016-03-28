/*
	Author Kardå Chalak & Andreas Hahr
	{kardac, "andreas adress"}@kth.se
	
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
	
	/*Andreas! I interpret the grid size as the length of the side.
	  Otherwise we have to take the sqrt of it to get the sides when
	  malloc'ing, change it if you see it fitting. :) */
	gridSize = (argc > 1)? atoi(argv[1]) + 2 : GRIDSIZE;
	numIters = (argc > 2)? atoi(argv[2]) * 0.5 : NUMITERS;
	
	printf("gridSize: %d\n", gridSize);
	printf("numIters: %d\n", numIters);
	
	/*malloc the columns*/
	float **grid = (float **) malloc(gridSize * sizeof(float*));
	float **new = (float **) malloc(gridSize * sizeof(float*));
	/*malloc the rows*/
	for(i = 0; i < gridSize; i++){
		grid[i] = (float *) malloc(gridSize * sizeof(float));
		new[i] = (float *) malloc(gridSize * sizeof(float));
	}
	/*Initiation of grid boundary */
	/*Hörnen initieras 2 gånger. Går säkert att göra bättre men just nu känns det
	  inte värt*/
	for(i = 0; i < gridSize; i++){
		grid[0][i] = 1;				// upper boundary
		new[0][i] =1;
		grid[gridSize-1][i] = 1;	// lower boundary
		new[gridSize-1][i] = 1;
		grid[i][0] = 1; 			// left  boundary
		new[i][0] = 1;
		grid[i][gridSize-1] =1;		// right boundary
		new[i][gridSize-1] = 1;
	}
	/*Instantiation of grids */
	for(i = 1; i < gridSize-1; i++){
		for(j = 1; j < gridSize-1; j++){
			grid[i][j] = 0;
			//new[i][j] = 0; // behövs nog inte
		}	
	}
	
	/*Jacobi iterations*/
	/*Körde * 0.625 på andra iterationen istället för 2st * 0.25
	  fick annat resultat, måste gjort fel, låter vara såhär så länge*/
	for(k = 0; k < numIters; k++){
	
		for(i = 1; i < gridSize-1; i++){
			for(j = 1; j < gridSize-1; j++){
				new[i][j] = (grid[i-1][j] +
							 grid[i+1][j] +
							 grid[i][j-1] +
							 grid[i][j+1])* 0.25 ;
			}	
		}
		
		for(i = 1; i < gridSize-1; i++){
			for(j = 1; j < gridSize-1; j++){
				grid[i][j] = (new[i-1][j] +
							  new[i+1][j] +
							  new[i][j-1] +
							  new[i][j+1])* 0.25;
			}
		}
	}
	/*Prints the grid with boundary*/
	for(i = 0; i < gridSize; i++){
			printf("\n");
		for(j = 0; j < gridSize; j++){
			printf("%f  ", grid[i][j]);
		}	
	}
	printf("\n");
}
