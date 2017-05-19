/*
A sequential PDE-solver for Laplace's equation using Jacobi iterative and
multigrid methods

Authors:  Andreas Hahr    Kardå Chalak
          hahr@kth.se     kardac@kth.se

arg1: gridSize, the grid size for the coarsest grid, not including boundaries
the true grid size will be 8xgridSize + 9, default = 24

arg2: numIters, nr. Jacobi iterations for coarsest grid, default = 10

arg3: finIter, nr. Jacobi iterations for finer grids, default = 4

Grids:  grid (coarsest),
        grid2,
        grid3,
        grid4 (finest)
*/

#include <stdlib.h>
#include <stdio.h>
#include <omp.h> // for using omp_get_wtime()

void jacobi(double** grid, double** new, int gridSize, int iter);
void restriction(double** fine, double** coarse, int coarseSize);
void interpolation(double** coarse, double** fine, int coarseSize, int fineSize);
void gridprint(double** grid, int gridSize);

int main(int argc, char *argv[]) {
  int gridSize2, gridSize3, gridSize4;
  double temp, maxdiff, timestamp = 0;
  int gridSize = 24, numIters = 10, finIter = 4;

  switch (argc) {
    case 4: finIter = atoi(argv[3]);
    case 3: numIters = atoi(argv[2]);
    case 2: gridSize = atoi(argv[1]);
  }

  gridSize2 = 2*gridSize+1;
  gridSize3 = 2*gridSize2+1;
  gridSize4 = 2*gridSize3+1;

  printf("gridSize: %d\n", gridSize);
  printf("numIters: %d\n", numIters);
  printf("finIter: %d\n", finIter);

  /* allocating memory for the grids */
  /* rows */
  double **grid = (double **) malloc((gridSize+2) * sizeof(double*));
  double **new = (double **) malloc((gridSize+2) * sizeof(double*));

  double **grid2 = (double **) malloc((gridSize2+2) * sizeof(double*));
  double **new2 = (double **) malloc((gridSize2+2) * sizeof(double*));

  double **grid3 = (double **) malloc((gridSize3+2) * sizeof(double*));
  double **new3 = (double **) malloc((gridSize3+2) * sizeof(double*));

  double **grid4 = (double **) malloc((gridSize4+2) * sizeof(double*));
  double **new4 = (double **) malloc((gridSize4+2) * sizeof(double*));

  /* columns */
  for(int i = 0; i < gridSize+2; i++) {
    grid[i] = (double *) malloc((gridSize+2) * sizeof(double));
    new[i] = (double *) malloc((gridSize+2) * sizeof(double));
  }
  for(int i = 0; i < gridSize2+2; i++) {
    grid2[i] = (double *) malloc((gridSize2+2) * sizeof(double));
    new2[i] = (double *) malloc((gridSize2+2) * sizeof(double));
  }
  for(int i = 0; i < gridSize3+2; i++) {
    grid3[i] = (double *) malloc((gridSize3+2) * sizeof(double));
    new3[i] = (double *) malloc((gridSize3+2) * sizeof(double));
  }
  for(int i = 0; i < gridSize4+2; i++) {
    grid4[i] = (double *) malloc((gridSize4+2) * sizeof(double));
    new4[i] = (double *) malloc((gridSize4+2) * sizeof(double));
  }

  /* set grid boundary points */
  for(int i = 0; i < gridSize+2; i++) {
    grid[0][i] = 1;				     // upper boundary
    new[0][i] = 1;
    grid[gridSize+1][i] = 1;	 // lower boundary
    new[gridSize+1][i] = 1;
    grid[i][0] = 1; 			     // left  boundary
    new[i][0] = 1;
    grid[i][gridSize+1] = 1;	 // right boundary
    new[i][gridSize+1] = 1;
  }

  /* set grid2 boundary points */
  for(int i = 0; i < gridSize2+2; i++) {
    grid2[0][i] = 1;            // upper boundary
    new2[0][i] = 1;
    grid2[gridSize2+1][i] = 1;  // lower boundary
    new2[gridSize2+1][i] = 1;
    grid2[i][0] = 1; 			      // left  boundary
    new2[i][0] = 1;
    grid2[i][gridSize2+1] = 1;  // right boundary
    new2[i][gridSize2+1] = 1;
  }

  /* set grid3 boundary points */
  for(int i = 0; i < gridSize3+2; i++) {
    grid3[0][i] = 1;            // upper boundary
    new3[0][i] = 1;
    grid3[gridSize3+1][i] = 1;  // lower boundary
    new3[gridSize3+1][i] = 1;
    grid3[i][0] = 1; 			      // left  boundary
    new3[i][0] = 1;
    grid3[i][gridSize3+1] = 1;  // right boundary
    new3[i][gridSize3+1] = 1;
  }

  /* set grid4 boundary points */
  for(int i = 0; i < gridSize4+2; i++) {
    grid4[0][i] = 1;            // upper boundary
    new4[0][i] = 1;
    grid4[gridSize4+1][i] = 1;  // lower boundary
    new4[gridSize4+1][i] = 1;
    grid4[i][0] = 1; 			      // left  boundary
    new4[i][0] = 1;
    grid4[i][gridSize4+1] = 1;  // right boundary
    new4[i][gridSize4+1] = 1;
  }

  /* set grid4 interior points */
  for(int i = 1; i < gridSize4+1; i++) {
    for(int j = 1; j < gridSize4+1; j++) {
      grid4[i][j] = 0;
    }
  }

  timestamp -= omp_get_wtime(); // time

  /* start V-cycle */
  jacobi(grid4, new4, gridSize4, finIter);
  restriction(grid4, grid3, gridSize3);

  jacobi(grid3, new3, gridSize3, finIter);
  restriction(grid3, grid2, gridSize2);

  jacobi(grid2, new2, gridSize2, finIter);
  restriction(grid2, grid, gridSize);

  jacobi(grid, new, gridSize, numIters);

  interpolation(grid, grid2, gridSize, gridSize2);
  jacobi(grid2, new2, gridSize2, finIter);

  interpolation(grid2, grid3, gridSize2, gridSize3);
  jacobi(grid3, new3, gridSize3, finIter);

  interpolation(grid3, grid4, gridSize3, gridSize4);
  jacobi(grid4, new4, gridSize4, finIter);
  /* V-cycle ended */

  timestamp += omp_get_wtime(); // time

  gridprint(grid4, gridSize4); // print grid to filedata.out

  /* computation of the maximum difference */
  /*maxdiff = 0;
  for(int i = 1; i < gridSize4-1; i++) {
    for(int j = 1; j < gridSize4-1; j++) {
      temp = grid4[i][j] - new4[i][j];
      if(temp < 0) temp = -temp;
      if(temp > maxdiff) maxdiff = temp;
    }
  }
  printf("Maximum difference: %f\n", maxdiff);*/

  printf("Maximum final error: %f\n", 1.0-grid4[(gridSize4+1)/2][(gridSize4+1)/2]);
  printf("Jacobi iterations runtime: %f seconds\n", timestamp);

  exit(0);
}


void jacobi(double** grid, double** new, int gridSize, int iter) {
  for(int k = 0; k < iter; k++) {
    for(int i = 1; i < gridSize+1; i++) {
      for(int j = 1; j < gridSize+1; j++) {
        new[i][j] = (grid[i-1][j] // value from left
        + grid[i+1][j]            // value from right
        + grid[i][j-1]            // value from above
        + grid[i][j+1]) * 0.25;   // value from below
      }
    }
    for(int i = 1; i < gridSize+1; i++) {
      for(int j = 1; j < gridSize+1; j++) {
        grid[i][j] = (new[i-1][j] // value from left
        + new[i+1][j]             // value from right
        + new[i][j-1]             // value from above
        + new[i][j+1]) * 0.25;    // value from below
      }
    }
  }
}

void restriction(double** fine, double** coarse, int coarseSize) {
  for(int i = 1; i < coarseSize+1 ; i++) {
    int m = i << 1;
    for(int j = 1; j < coarseSize+1; j++) {
      int n = j << 1;
      coarse[i][j] = fine[m][n] * 0.5 +
      (fine[m-1][n] + fine[m+1][n] + fine[m][n-1] + fine[m][n+1]) * 0.125;
    }
  }
}

void interpolation(double** coarse, double** fine, int coarseSize, int fineSize) {
  /* place coarse grid points on finer grid */
  for(int i = 1; i < coarseSize+1; i++) {
    int m = i << 1;
    for(int j = 1; j < coarseSize+1; j++) {
      int n = j << 1;
      fine[m][n] = coarse[i][j];
    }
  }
  /* set fine grid points in columns corresponding to coarser grid */
  for(int j = 2; j < fineSize+1; j+=2) {
    for(int i = 1; i < fineSize+1; i+=2) {
      fine[i][j] = (fine[i-1][j] + fine[i+1][j]) * 0.5;
    }
  }
  /* set rest of the points */
  for(int i = 1; i < fineSize+1; i++) {
    for(int j = 1; j < fineSize+1; j+=2) {
      fine[i][j] = (fine[i][j-1] + fine[i][j+1]) * 0.5;
    }
  }
}

void gridprint(double** grid, int gridSize) {
  FILE *fp;
  fp = fopen("filedata.out", "w+");
  for(int i = 0; i < gridSize+2; i++) {
    for(int j = 0; j < gridSize+2; j++) {
      fprintf(fp, "%f  ", grid[i][j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}
