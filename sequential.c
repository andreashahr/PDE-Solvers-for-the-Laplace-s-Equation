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
	
	int gridSize;
	int numIters;
	
	gridSize = (argc > 1)? atoi(argv[1]) : GRIDSIZE;
	numIters = (argc > 2)? atoi(argv[2]) : NUMITERS;
	
	printf("gridSize: %d\n", gridSize);
	printf("numIters: %d\n", numIters);
}
