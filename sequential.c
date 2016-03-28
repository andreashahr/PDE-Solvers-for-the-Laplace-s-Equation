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
	
	int gridSize, numIters, i;
	
	/*Andreas I interpret the grid size as the length of the side.
	  Otherwise we have to take the sqrt of it to get the sides when
	  malloc'ing, change it if you see it fitting. :) */
	gridSize = (argc > 1)? atoi(argv[1]) : GRIDSIZE;
	numIters = (argc > 2)? atoi(argv[2]) : NUMITERS;
	
	printf("gridSize: %d\n", gridSize);
	printf("numIters: %d\n", numIters);
	
	/*malloc the columns*/
	int **grid = (int **) malloc(gridSize+2 * sizeof(int*));
	int **new = (int **) malloc(gridSize+2 * sizeof(int*));
	/*malloc the rows*/
	for(i = 0; i < gridSize+2; i++){
		grid[i] = (int *) malloc(gridSize+2 * sizeof(int));
		new[i] = (int *) malloc(gridSize+2 * sizeof(int));
	}
}
