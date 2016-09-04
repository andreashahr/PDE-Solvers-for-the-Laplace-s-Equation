/*
	Author Kardå Chalak & Andreas Hahr
	kardac@kth.se
	andreas epost
	
	A sequential Jacobi Multigrid iteration program
	Grids are: gridh (finest),
			   grid2h,
			   grid4h,
			   grid8h (coarsest).
*/

/*

Example of a gridSize of 8 (+ boundary)  and its multigrids

									B  B  B  B  B  B  B  B  
								 B  *  *  *  *  *  *  *  *  B
								 B  *  *  *  *  *  *  *  *  B
								 B  *  *  *  *  *  *  *  *  B
								 B  *  *  *  *  *  *  *  *  B
								 B  *  *  *  *  *  *  *  *  B
								 B  *  *  *  *  *  *  *  *  B
								 B  *  *  *  *  *  *  *  *  B
								 B  *  *  *  *  *  *  *  *  B
									B  B  B  B  B  B  B  B  
 
 
    B  B  B  B  B  B  B  B  		 B    B    B    B
 B  0  *  0  *  0  *  0  *  B	   B *    *    *    *    B
 B  *  *  *  *  *  *  *  *  B 		
 B  0  *  0  *  0  *  0  *  B	   B *    *    *    *    B  	
 B  *  *  *  *  *  *  *  *  B		
 B  0  *  0  *  0  *  0  *  B      B *    *    *    *    B 
 B  *  *  *  *  *  *  *  *  B		
 B  0  *  0  *  0  *  0  *  B      B *    *    *    *    B 
 B  *  *  *  *  *  *  *  *  B		
    B  B  B  B  B  B  B  B           B    B    B    B
	
 
    B  B  B  B  	     			   B    B
 B  0  *  0  *  B	  				B  *    *  B
 B  *  *  *  *  B
 B  0  *  0  *  B 				    B  *    *  B
 B  *  *  *  *  B
    B  B  B  B  			           B    B
 

    B  B  							B
 B  0  *  B     				 B  *     B
 B  *  *  B
    B  B  				            B

 
 Example of a grid size of 9 (+ boundary) and its multigrids:
 
										B  B  B  B  B  B  B  B  B
									 B  *  *  *  *  *  *  *  *  *  B
									 B  *  *  *  *  *  *  *  *  *  B
									 B  *  *  *  *  *  *  *  *  *  B
									 B  *  *  *  *  *  *  *  *  *  B
									 B  *  *  *  *  *  *  *  *  *  B
									 B  *  *  *  *  *  *  *  *  *  B
									 B  *  *  *  *  *  *  *  *  *  B
									 B  *  *  *  *  *  *  *  *  *  B
									 B  *  *  *  *  *  *  *  *  *  B
										B  B  B  B  B  B  B  B  B
 
 
    B  B  B  B  B  B  B  B  B		     B     B     B     B     B	
 B  0  *  0  *  0  *  0  *  0  B      B  *     *     *     *     *  B 
 B  *  *  *  *  *  *  *  *  *  B
 B  0  *  0  *  0  *  0  *  0  B      B  *     *     *     *     *  B
 B  *  *  *  *  *  *  *  *  *  B
 B  0  *  0  *  0  *  0  *  0  B      B  *     *     *     *     *  B
 B  *  *  *  *  *  *  *  *  *  B
 B  0  *  0  *  0  *  0  *  0  B      B  *     *     *     *     *  B
 B  *  *  *  *  *  *  *  *  *  B
 B  0  *  0  *  0  *  0  *  0  B      B  *     *     *     *     *  B 
    B  B  B  B  B  B  B  B  B            B     B     B     B     B
	

   B  B  B  B  B          					  B     B     B
B  0  *  0  *  0  B     				   B  *     *     *  B
B  *  *  *  *  *  B      
B  0  *  0  *  0  B     				   B  *     *     *  B
B  *  *  *  *  *  B 
B  0  *  0  *  0  B    					   B  *     *     *  B
   B  B  B  B  B          					  B     B     B
 
 
   B  B  B									B     B
 B 0  *  0  B						    B	*     *  B
 B *  *  *  B		  
 B 0  *  0  B		 				 	B	*     *  B
   B  B  B									B     B
   
We can see that we need different restriction and interpolation operators for even and uneven
grids. When the grid size is even the on step down courser grid is computed by gridsize/2 +1.
When the grid is uneven it is computed by gridsize/2 + 2. 

However when implementing the restriction operator method I have only taken the uneven case
into consideration. This is improvable. The uneven case is choose as it is more
trivial to cover the same spatial space with a courser grid. 
  
 Big risk for false sharing in cashe as the program uses a lot of matrix calculation
 hur ställer du dig till detta andreas, är det en risk, bör vi försöka lösa det redan nu
 eller först göra prestanda tester och se om det behövs?
*/
 #include <stdlib.h>
 #include <stdio.h>
 
 #define GRIDSIZE 10
 #define NUMITERS 10
 #define FINERITER 2
 
 void jacobi(double** grid, double** new, int gridSize, int iter);
 void restriction(double** fine, double** coarse, int fineSize, int coarseSize);
 void interpolation(double** coarse, double** fine, int coarseSize, int fineSize);
 
 void gridprint(double** grid, int gridSize);
 
 int main(int argc, char *argv[]){
	/*gridSize, the grid size for the finest grid, not including boundaries
				the grid size has to be at least 8 or bigger
				
	  numIters, the number of iterations to use for the coarsest grid
				the finer grids has 4 iterations each (specification from Prof. Dr. Vladimir Vlassov). 
	*/
	 
	printf("First print\n"); 
	int gridSizeH, gridSize2H, gridSize4H, gridSize8H, numIters, i, j, k;
	double temp, maxdiff;
	 
	gridSize8H = (argc > 1)? atoi(argv[1]) : GRIDSIZE;
	gridSize4H = gridSize8H * 2 + 1;
	gridSize2H = gridSize4H * 2 + 1;
	gridSizeH = gridSize2H * 2 + 1;
	
	numIters = (argc > 2)? atoi(argv[2]) : NUMITERS;
	 
	printf("gridSizeH: %d\n", gridSizeH);
	printf("gridSize2H: %d\n", gridSize2H);
	printf("gridSize4H: %d\n", gridSize4H);
	printf("gridSize8H: %d\n", gridSize8H);
	printf("numIters: %d\n", numIters);
	
	/*Allocating memory for the grids*/
	
	/*malloc the rows*/
	double **gridH = (double **) malloc((gridSizeH + 2 )* sizeof(double*));
	double **newH = (double **) malloc((gridSizeH + 2) * sizeof(double*));
	
	double **grid2H = (double **) malloc((gridSize2H + 2) * sizeof(double*));
	double **new2H = (double **) malloc((gridSize2H + 2) * sizeof(double*));

	double **grid4H = (double **) malloc((gridSize4H + 2) * sizeof(double*));
	double **new4H = (double **) malloc((gridSize4H + 2) * sizeof(double*));
	
	double **grid8H = (double **) malloc((gridSize8H + 2 ) * sizeof(double*));
	double **new8H = (double **) malloc((gridSize8H + 2) * sizeof(double*));
	
	
	/*malloc the columns*/
	for(i = 0; i < gridSizeH + 2; i++){
		gridH[i] = (double *) malloc((gridSizeH + 2) * sizeof(double));
		newH[i] = (double *) malloc((gridSizeH + 2) * sizeof(double));
	}
	for(i = 0; i < gridSize2H + 2; i++){
		grid2H[i] = (double *) malloc((gridSize2H + 2) * sizeof(double));
		new2H[i] = (double *) malloc((gridSize2H + 2) * sizeof(double));
	}
	for(i = 0; i < gridSize4H + 2; i++){
		grid4H[i] = (double *) malloc((gridSize4H + 2) * sizeof(double));
		new4H[i] = (double *) malloc((gridSize4H + 2) * sizeof(double));
	}
	for(i = 0; i < gridSize8H + 2; i++){
		grid8H[i] = (double *) malloc((gridSize8H + 2) * sizeof(double));
		new8H[i] = (double *) malloc((gridSize8H + 2 ) * sizeof(double));
	}
	
	/*Initiation of gridH boundary */
	for(i = 0; i < gridSizeH + 2; i++){
		gridH[0][i] = 1;				// upper boundary
		newH[0][i] = 1;
		gridH[gridSizeH+1][i] = 1;	// lower boundary
		newH[gridSizeH+1][i] = 1;
		gridH[i][0] = 1; 			// left  boundary
		newH[i][0] = 1;
		gridH[i][gridSizeH+1] = 1;	// right boundary
		newH[i][gridSizeH+1] = 1;
	}
	
	/*Initiation of grid2H boundary */
	for(i = 0; i < gridSize2H + 2; i++){
		grid2H[0][i] = 1;				// upper boundary
		new2H[0][i] = 1;
		grid2H[gridSize2H+1][i] = 1;	// lower boundary
		new2H[gridSize2H+1][i] = 1;
		grid2H[i][0] = 1; 			// left  boundary
		new2H[i][0] = 1;
		grid2H[i][gridSize2H+1] = 1;	// right boundary
		new2H[i][gridSize2H+1] = 1;
	} 
	
	/*Initiation of grid4H boundary */
	for(i = 0; i < gridSize4H + 2; i++){
		grid4H[0][i] = 1;				// upper boundary
		new4H[0][i] = 1;
		grid4H[gridSize4H+1][i] = 1;	// lower boundary
		new4H[gridSize4H+1][i] = 1;
		grid4H[i][0] = 1; 			// left  boundary
		new4H[i][0] = 1;
		grid4H[i][gridSize4H+1] = 1;	// right boundary
		new4H[i][gridSize4H+1] = 1;
	}
	
		/*Initiation of grid8H boundary */
	for(i = 0; i < gridSize8H + 2; i++){
		grid8H[0][i] = 1;				// upper boundary
		new8H[0][i] = 1;
		grid8H[gridSize8H+1][i] = 1;	// lower boundary
		new8H[gridSize8H+1][i] = 1;
		grid8H[i][0] = 1; 			// left  boundary
		new8H[i][0] = 1;
		grid8H[i][gridSize8H+1] = 1;	// right boundary
		new8H[i][gridSize8H+1] = 1;
	}
	
	/*Instantiation of gridH */
	for(i = 1; i < gridSizeH+1; i++){
		for(j = 1; j < gridSizeH+1; j++){
			gridH[i][j] = 0;
		}	
	}
	
	/*Jacobi iterations on gridH
	  FINERITER * 2 = number of iterations
	*/
	
	printf("starting V iterations");
	gridprint(gridH, gridSizeH);
	jacobi(gridH, newH, gridSizeH, FINERITER);
	gridprint(gridH, gridSizeH);

	printf("before first restriction");	
	restriction(gridH, grid2H, gridSizeH, gridSize2H);
	gridprint(grid2H, gridSize2H);
	jacobi(grid2H, new2H, gridSize2H, FINERITER);
	gridprint(grid2H, gridSize2H);
	
	printf("berfore secound restriciton\n");
	restriction(grid2H, grid4H, gridSize2H, gridSize4H);
	gridprint(grid4H, gridSize4H);
	jacobi(grid4H, new4H, gridSize4H, FINERITER);
	gridprint(grid4H, gridSize4H);

	printf("berfore thrird restriction\n");	
	restriction(grid4H, grid8H, gridSize4H, gridSize8H);
	printf("sista restriciton\n");
	gridprint(grid8H, gridSize8H);
	printf("innan jacobi\n");
	jacobi(grid8H, new8H, gridSize8H, numIters);
	printf("efter jacobi\n");
	gridprint(grid8H, gridSize8H);
	
	printf("innan interpolation\n");
	interpolation(grid8H, grid4H, gridSize8H, gridSize4H);
	printf("interpolation lyckades \n");
	gridprint(grid4H, gridSize4H);
	jacobi(grid4H, new4H, gridSize4H, FINERITER);
	gridprint(grid4H, gridSize4H);
	
	interpolation(grid4H, grid2H, gridSize4H, gridSize2H);
	gridprint(grid2H, gridSize2H);
	jacobi(grid2H, new2H, gridSize2H, FINERITER);
	gridprint(grid2H, gridSize2H);
	
	interpolation(grid2H, gridH, gridSize2H, gridSizeH);
	gridprint(gridH, gridSizeH);
	jacobi(gridH, newH, gridSizeH, FINERITER);
	gridprint(gridH, gridSizeH);
	
	
	free(gridH);
	free(grid2H);
	free(grid4H);
	free(grid8H);
	free(newH);
	free(new2H);
	free(new4H);
	free(new8H);
 }
 	/*Vi bör eventuellt använda shift istället för * 0.25*/
 void jacobi(double** grid, double** new, int gridSize, int iter){
 
	int i, j, k;
	
	for(k = 0; k < iter; k++) {
	
		for(i = 1; i < gridSize+1; i++) {
			for(j = 1; j < gridSize+1; j++) {
				new[i][j] = (grid[i-1][j] +
							 grid[i+1][j] +
							 grid[i][j-1] +
							 grid[i][j+1]) * 0.25;
			}
		}
		for(i = 1; i < gridSize+1; i++) {
			for(j = 1; j < gridSize+1; j++) {
				grid[i][j] = (new[i-1][j] +
							  new[i+1][j] +
							  new[i][j-1] +
							  new[i][j+1]) * 0.25;
			}
		}
	}
 }
 
 /*
	Restriction operator uses the following transformation matrices:
	(B stands for boundary value and its corresponding value is 0)
	
 */
 void restriction(double** fine, double** coarse, int fineSize, int coarseSize){
	
	/*
	Restriction matrix 
	
	0		1/8		0
	
	1/8		1/2		1/8
	
	0		1/8		0
	
	*/
	
	int i, j, m, n;
	
	for(i = 1; i < coarseSize+1; i++){
		m = i * 2;
		
		for(j = 1; j < coarseSize+1; j++){
			n = j * 2;
			coarse[i][j] = fine[m][n] * 0.5 + (fine[m-1][n] + fine[m][n-1] + 
				                         fine[m][n+1] + fine[m+1][n]) * 0.125;
		}
	}
	
 }
void interpolation(double** coarse, double** fine, int coarseSize, int fineSize){
 
	/*
		Uses the interpolation matrix
		
		1/4		1/2		1/4
		
		1/2		1		1/2
		
		1/4		1/2		1/4
	*/
	
	int i, j, n, m;	
	// The centre part (the 1) of the interpolation matrix
	for(i = 1; i < coarseSize+2; i++){
		// adapting to the finer matrix positions
		m = i * 2;	
		
		for(j = 1; j < coarseSize+2; j++){
			// adapting to the coarser matrix positions
			n = j * 2;
			
			fine[m][n] = coarse[i][j];
		}	
	}
	/*The surrounding elements in the matrix */
	// initiliazise columns that has coresponding coarse values
	for(j = 2; j <fineSize+2; j += 2){
		for(i =1; i <fineSize+2; i += 2){
			fine[i][j] = (fine[i-1][j] + fine[i+1][j]) * 0.5;
		}
	}
	// initialize the rest of the columns
	for(j =1; j<fineSize+2; j+=2){
		for(i=1; i<fineSize+2; i++){
			fine[i][j] = (fine[i][j-1] + fine[i][j+1]) *0.5;
		}
	}
 }
 
 void gridprint(double** grid, int gridSize){
 
	int i, j;
	
	for(i = 0; i < gridSize+2; i++){
		printf("\n");
		for(j = 0; j < gridSize+2; j++){
			printf("%f  ", grid[i][j]);
		}
	}
	printf("\n");
 }
