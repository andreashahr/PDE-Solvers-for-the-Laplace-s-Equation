/*
A parallel PDE-solver for Laplace's equation using Jacobi iterative method

Authors:  Andreas Hahr    Kardå Chalak
          hahr@kth.se     kardac@kth.se

arg1: gridSize (int), the grid side length with boundary points excluded,
default = 10

arg2: numIters (int), nr. Jacobi iterations, default = 10 (20)

arg3: numProc (int), nr. of threads to use, default = omp_get_max_threads()
*/

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

void dissBarrier(int* counters);
void gridprint(double** grid, int gridSize);

int numProc;

int main(int argc, char *argv[]) {
  int gridSize = 10, numIters = 10;
  double temp, maxdiff, timestamp = 0;
  numProc = omp_get_max_threads();

  switch (argc) {
    case 4: numProc = atoi(argv[3]);
    case 3: numIters = atoi(argv[2]);
    case 2: gridSize = atoi(argv[1]);
  }

  /* print option settings */
  printf("gridSize: %d\n", gridSize);
  printf("numIters: %d\n", numIters);
  printf("numThreads: %d\n", numProc);

  gridSize += 2;

  /* allocate memory for barrier flags */
  int* barFlags = (int *) calloc(numProc, sizeof(int));

  /* allocate memory for grids */
  /* columns */
  double **grid = (double **) malloc(gridSize * sizeof(double*));
  double **new = (double **) malloc(gridSize * sizeof(double*));
  /* rows */
  for(int i = 0; i < gridSize; i++) {
    grid[i] = (double *) malloc(gridSize * sizeof(double));
    new[i] = (double *) malloc(gridSize * sizeof(double));
  }

  /* set grid boundary */
  for(int i = 0; i < gridSize; i++) {
    grid[0][i] = 1;				     // upper boundary
    new[0][i] = 1;
    grid[gridSize-1][i] = 1;	 // lower boundary
    new[gridSize-1][i] = 1;
    grid[i][0] = 1;			       // left  boundary
    new[i][0] = 1;
    grid[i][gridSize-1] = 1;	 // right boundary
    new[i][gridSize-1] = 1;
  }

  /* set interior points */
  for(int i = 1; i < gridSize-1; i++) {
    for(int j = 1; j < gridSize-1; j++) {
      grid[i][j] = 0;
    }
  }

  timestamp -= omp_get_wtime(); // Jacobi iterations start timestamp

  #pragma omp parallel num_threads(numProc) // fork
  {
  /* Jacobi iterations */
  for(int k = 0; k < numIters; k++) {
    #pragma omp for nowait schedule(static)
    for(int i = 1; i < gridSize-1; i++) {
      for(int j = 1; j < gridSize-1; j++) {
        new[i][j] = (grid[i-1][j]   // value from left
        + grid[i+1][j]              // value from right
        + grid[i][j-1]              // value from below
        + grid[i][j+1]) * 0.25;     // value from above
      }
    }
    dissBarrier(barFlags);
    #pragma omp for nowait schedule(static)
    for(int i = 1; i < gridSize-1; i++) {
      for(int j = 1; j < gridSize-1; j++) {
        grid[i][j] = (new[i-1][j]   // value from left
        + new[i+1][j]               // value from right
        + new[i][j-1]               // value from below
        + new[i][j+1]) * 0.25;      // value from above
      }
    }
    dissBarrier(barFlags);
  }
  } // join threads

  timestamp += omp_get_wtime(); // Jacobi iterations end timestamp

  gridprint(grid, gridSize); // print grid to filedata.out

  /* computation of the maximum difference */
  /*maxdiff = 0;
  for(int i = 1; i < gridSize-1; i++) {
    for(int j = 1; j < gridSize-1; j++) {
      temp = grid[i][j] - new[i][j];
      if(temp < 0) temp = -temp;
      if(temp > maxdiff) maxdiff = temp;
    }
  }
  printf("Maximum difference: %f\n", maxdiff);*/

  printf("Jacobi iterations runtime: %f seconds\n", timestamp);
  printf("Maximum final error: %f\n", 1.0-grid[(gridSize+1)/2][(gridSize+1)/2]);

  exit(0);
}

void dissBarrier(int* barFlags) {
  int me = omp_get_thread_num();
  for(int s = 1; s < numProc; s <<= 1) {
    barFlags[me] += 1;
    while(barFlags[me] > barFlags[(me+s)%numProc]) {}; // wait
  }
}

void gridprint(double** grid, int gridSize) {
  FILE *fp;
  fp = fopen("filedata.out", "w+");
  for(int i = 0; i < gridSize; i++) {
    for(int j = 0; j < gridSize; j++) {
      fprintf(fp, "%f  ", grid[i][j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}
