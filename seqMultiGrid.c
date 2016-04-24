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
  
	
*/
 #include <stdlib.h>
 #include <stdio.h>
 
 #define GRIDSIZE 10
 #define NUMITERS 10
 
 int main(int argc, char *argv[]){
	/*gridSize, the grid size for the finest grid, not including boundaries
				the grid size has to be at least 8 or bigger
				
	  numIters, the number of iterations to use for the coarsest grid
				the finer grids has 4 iterations each (specification from Prof. Dr. Vladimir Vlassov). 
	*/
	 
	int gridSizeH, gridSize2H, gridSize4H, gridSize8H, numIters, i, j, k;
	double temp, maxdiff;
	 
	gridSizeH = (argc > 1)? atoi(argv[1]) + 2 : GRIDSIZE;
	if(gridSizeH % 2 ){
		// uneven
		gridSize2H = gridSizeH * 0.5 + 2;
	}else{
		//even
		gridSize2H = gridSizeH *0.5 + 1;
	}
	if(gridSize2H % 2 ){
		// uneven
		gridSize4H = gridSize2H * 0.5 + 2;
	}else{
		//even
		gridSize4H = gridSize2H *0.5 + 1;
	}
	if(gridSizeH % 2 ){
		// uneven
		gridSize8H = gridSize4H * 0.5 + 2;
	}else{
		//even
		gridSize8H = gridSize4H *0.5 + 1;
	}
	
	numIters = (argc > 2)? atoi(argv[2]) : NUMITERS;
	 
	printf("gridSizeH: %d\n", gridSizeH);
	printf("gridSize2H: %d\n", gridSize2H);
	printf("gridSize4H: %d\n", gridSize4H);
	printf("gridSize8H: %d\n", gridSize8H);
	printf("numIters: %d\n", numIters);
	
	/*Allocating memory for the grids*/
	
	/*malloc the columns*/
	double **gridH = (double **) malloc(gridSizeH * sizeof(double*));
	double **newH = (double **) malloc(gridSizeH * sizeof(double*));
	
	double **grid2H = (double **) malloc(gridSize2H * sizeof(double*));
	double **new2H = (double **) malloc(gridSize2H * sizeof(double*));

	double **grid4H = (double **) malloc(gridSize4H * sizeof(double*));
	double **new4H = (double **) malloc(gridSize4H * sizeof(double*));
	
	double **grid8H = (double **) malloc(gridSize8H * sizeof(double*));
	double **new8H = (double **) malloc(gridSize8H * sizeof(double*));
	
	
	/*malloc the rows*/
	for(i = 0; i < gridSizeH; i++){
		gridH[i] = (double *) malloc(gridSizeH * sizeof(double));
		newH[i] = (double *) malloc(gridSizeH * sizeof(double));
	}
	for(i = 0; i < gridSize2H; i++){
		grid2H[i] = (double *) malloc(gridSize2H * sizeof(double));
		new2H[i] = (double *) malloc(gridSize2H * sizeof(double));
	}
	for(i = 0; i < gridSize4H; i++){
		grid4H[i] = (double *) malloc(gridSize4H * sizeof(double));
		new4H[i] = (double *) malloc(gridSize4H * sizeof(double));
	}
	for(i = 0; i < gridSize8H; i++){
		grid8H[i] = (double *) malloc(gridSize8H * sizeof(double));
		new8H[i] = (double *) malloc(gridSize8H * sizeof(double));
	}
	
	/*Initiation of gridH boundary */
	for(i = 0; i < gridSizeH; i++){
		gridH[0][i] = 1;				// upper boundary
		newH[0][i] = 1;
		gridH[gridSizeH-1][i] = 1;	// lower boundary
		newH[gridSizeH-1][i] = 1;
		gridH[i][0] = 1; 			// left  boundary
		newH[i][0] = 1;
		gridH[i][gridSizeH-1] = 1;	// right boundary
		newH[i][gridSizeH-1] = 1;
	}
	
	/*Initiation of grid2H boundary */
	for(i = 0; i < gridSize2H; i++){
		grid2H[0][i] = 1;				// upper boundary
		new2H[0][i] = 1;
		grid2H[gridSize2H-1][i] = 1;	// lower boundary
		new2H[gridSize2H-1][i] = 1;
		grid2H[i][0] = 1; 			// left  boundary
		new2H[i][0] = 1;
		grid2H[i][gridSize2H-1] = 1;	// right boundary
		new2H[i][gridSize2H-1] = 1;
	} 
	
	/*Initiation of grid4H boundary */
	for(i = 0; i < gridSize4H; i++){
		grid4H[0][i] = 1;				// upper boundary
		new4H[0][i] = 1;
		grid4H[gridSize4H-1][i] = 1;	// lower boundary
		new4H[gridSize4H-1][i] = 1;
		grid4H[i][0] = 1; 			// left  boundary
		new4H[i][0] = 1;
		grid4H[i][gridSize4H-1] = 1;	// right boundary
		new4H[i][gridSize4H-1] = 1;
	}
	
		/*Initiation of grid8H boundary */
	for(i = 0; i < gridSize8H; i++){
		grid8H[0][i] = 1;				// upper boundary
		new8H[0][i] = 1;
		grid8H[gridSize8H-1][i] = 1;	// lower boundary
		new8H[gridSize8H-1][i] = 1;
		grid8H[i][0] = 1; 			// left  boundary
		new8H[i][0] = 1;
		grid8H[i][gridSize8H-1] = 1;	// right boundary
		new8H[i][gridSize8H-1] = 1;
	}
	
	/*Instantiation of gridH */
	for(i = 1; i < gridSizeH-1; i++){
		for(j = 1; j < gridSizeH-1; j++){
			gridH[i][j] = 0;
		}	
	}
	
	/*Instantiation of grid2H */
	for(i = 1; i < gridSize2H-1; i++){
		for(j = 1; j < gridSize2H-1; j++){
			grid2H[i][j] = 0;
		}	
	}
	
	/*Instantiation of grid4H */
	for(i = 1; i < gridSize4H-1; i++){
		for(j = 1; j < gridSize4H-1; j++){
			grid4H[i][j] = 0;
		}	
	}
	
	/*Instantiation of grid8H */
	for(i = 1; i < gridSize8H-1; i++){
		for(j = 1; j < gridSize8H-1; j++){
			grid8H[i][j] = 0;
		}	
	}
	
	/*test of multigrids instantiation */
	
	for(i = 0; i < gridSizeH; i++){
		printf("\n");
		for(j = 0; j < gridSizeH; j++){
			printf("%f  ", gridH[i][j]);
		}
	}
	printf("\n");
	
	for(i = 0; i < gridSize2H; i++){
		printf("\n");
		for(j = 0; j < gridSize2H; j++){
			printf("%f  ", grid2H[i][j]);
		}
	}
	printf("\n");

	for(i = 0; i < gridSize4H; i++){
		printf("\n");
		for(j = 0; j < gridSize4H; j++){
			printf("%f  ", grid4H[i][j]);
		}
	}
	printf("\n");
	
	for(i = 0; i < gridSize8H; i++){
		printf("\n");
		for(j = 0; j < gridSize8H; j++){
			printf("%f  ", grid8H[i][j]);
		}
	}
	printf("\n");
	
	free(gridH);
	free(grid2H);
	free(grid4H);
	free(grid8H);
	free(newH);
	free(new2H);
	free(new4H);
	free(new8H);
 }
 
 