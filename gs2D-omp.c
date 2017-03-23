
#ifdef _OPENMP
#include <omp.h>
#define RUN_TYPE "parallel"
#else
#define RUN_TYPE "serial"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "2Dgrid-verification.c"
#include "util.h"


void compute_gs_iterate(double*, double*, int);
void compute_gs_iterate_serial(double*, double*, int);
void compute_gs_iterate_parallel(double*, double*, int);
double* run_gs(int, double*, int);

int main(int argc, char **argv) {
    if(argc != 3) {
        printf("Usage: %s gridsize max_iterations\n", argv[0]);
        return EXIT_FAILURE;
    }
    int N = atoi(argv[1]);
    int max_iter = atoi(argv[2]);
    double *f = malloc(N * N * sizeof(double));
    double *u;
    int ind;
    double elapsed;
    timestamp_type t_start, t_end;
    for(ind = 0; ind < N * N; ind++) {
        f[ind] = 1;
    }
    get_timestamp(&t_start);
    u = run_gs(max_iter, f, N);
    get_timestamp(&t_end);
    elapsed = timestamp_diff_in_seconds(t_start, t_end);
    printf("Elapsed time was %.3f seconds on %dx%d gridpoints and  <=%d\
 iterations with a %s Gauss-Seidel method.\n", elapsed, N, N, max_iter, RUN_TYPE);
    free(u);
    free(f);
    return EXIT_SUCCESS;
}

double* run_gs(int max_iter, double* f, int N) {
    int iter;
    double rel_res, start_resid;
    double *u = calloc(N*N, sizeof(double));
    start_resid = compute_fd_residual(N, u, f);
    rel_res = 1;
    for(iter = 0; iter < max_iter && rel_res > 1e-4; iter++) {
        compute_gs_iterate(u, f, N);
        rel_res = compute_fd_residual(N, u, f) / start_resid;
    }
    printf("Final relative residual: %.5f\n", rel_res);
    return u;
}

void compute_gs_iterate(double* u, double* f, int N) {
    #ifdef _OPENMP
    compute_gs_iterate_parallel(u, f, N);
    #else
    compute_gs_iterate_serial(u, f, N);
    #endif
}

void compute_gs_iterate_serial(double* u, double* f, int N) {
    // Run a red-black GS method
    int i,j,color, ind;
    double h_sqr = 1.0 / ((N+1) * (N+1));
    for(color = 0; color < 2; color++) {
        //Handle the corners separately 
        if(color == 0) { //Red
            ind = 0;
            u[ind] = (h_sqr * f[ind] + u[ind+1] + u[ind+N]) / 4;
            ind = N*N - 1;
            u[ind] = (h_sqr * f[ind] + u[ind-1] + u[ind-N]) / 4;
            if(N%2 == 1) {
                //All corners red
                ind = N - 1;
                u[ind] = (h_sqr * f[ind] + u[ind - 1] + u[ind + N]) / 4;
                ind = N * (N-1);
                u[ind] = (h_sqr * f[ind] + u[ind + 1] + u[ind - N]) / 4;
            }
        } else { //Black
            if(N%2 == 0) {
                //Two black corners
                ind = N - 1;
                u[ind] = (h_sqr * f[ind] + u[ind - 1] + u[ind + N]) / 4;
                ind = N * (N-1);
                u[ind] = (h_sqr * f[ind] + u[ind + 1] + u[ind - N]) / 4;
            }
        }
        
        //Handle the top row
        for(ind = 2-color; ind < N - 1; ind+=2) {
            u[ind] = (h_sqr * f[ind] + u[ind - 1] + u[ind + 1]\
                    + u[ind + N]) / 4;
        }

        //Handle the bottom row
        for(ind = N*(N-1) + 2 - (color+N)%2; ind < N*N-1; ind+=2) {
            u[ind] = (h_sqr * f[ind] + u[ind - 1] + u[ind + 1] \
                        + u[ind - N]) / 4;
        }

        //Handle the interior rows
        for(i = 1; i < N-1; i++) {
            //Handle the edges
            if(i % 2 == color) {
                //Left edge
                ind = N*i;
                u[ind] = (h_sqr * f[ind] + u[ind+1] + u[ind+N] + u[ind-N]) / 4;
            }
            if((i + N - 1) % 2 == color) {
                //Right edge
                ind = N *  (i+1) - 1;
                u[ind] = (h_sqr * f[ind] + u[ind-1] + u[ind+N] + u[ind-N]) / 4;
            }
            for(j = 2 - (i+color)%2; j < N-1; j+=2) {
                ind = N*i + j;
                u[ind] = (h_sqr * f[ind] + u[ind-1] + u[ind+1] + u[ind+N]\
                            + u[ind-N]) / 4;
            }
        }
    }
}

void compute_gs_iterate_parallel(double* u, double* f, int N) {
    // Run a red-black GS method
    int i,j,color, ind;
    double h_sqr = 1.0 / ((N+1) * (N+1));
    for(color = 0; color < 2; color++) {
        #pragma omp parallel shared(u, f, N, color, h_sqr)\
            private(i,j,ind)
        {
            #pragma omp sections nowait
            {
                #pragma omp section
                {
                    //Handle the corners separately 
                    if(color == 0) { //Red
                        ind = 0;
                        u[ind] = (h_sqr * f[ind] + u[ind+1] + u[ind+N]) / 4;
                        ind = N*N - 1;
                        u[ind] = (h_sqr * f[ind] + u[ind-1] + u[ind-N]) / 4;
                        if(N%2 == 1) {
                            //All corners red
                            ind = N - 1;
                            u[ind] = (h_sqr * f[ind] + u[ind - 1] + u[ind + N])\
                                    / 4;
                            ind = N * (N-1);
                            u[ind] = (h_sqr * f[ind] + u[ind + 1] + u[ind - N])\
                                     / 4;
                        }
                    } else { //Black
                        if(N%2 == 0) {
                            //Two black corners
                            ind = N - 1;
                            u[ind] = (h_sqr * f[ind] + u[ind - 1] + u[ind + N])\
                                     / 4;
                            ind = N * (N-1);
                            u[ind] = (h_sqr * f[ind] + u[ind + 1] + u[ind - N])\
                                     / 4;
                        }
                    }
                }
                
                #pragma omp section
                {
                    //Handle the top row
                    for(ind = 2-color; ind < N - 1; ind+=2) {
                        u[ind] = (h_sqr * f[ind] + u[ind - 1] + u[ind + 1]\
                                + u[ind + N]) / 4;
                    }
                }

                #pragma omp section
                {
                    //Handle the bottom row
                    for(ind = N*(N-1) + 2 - (color+N)%2; ind < N*N-1; ind+=2) {
                        u[ind] = (h_sqr * f[ind] + u[ind - 1] + u[ind + 1] \
                                    + u[ind - N]) / 4;
                    }
                }
            }

            //Handle the interior rows
            #pragma omp for
            for(i = 1; i < N-1; i++) {
                //Handle the edges
                if(i % 2 == color) {
                    //Left edge
                    ind = N*i;
                    u[ind] = (h_sqr * f[ind] + u[ind+1] + u[ind+N] + u[ind-N])\
                            / 4;
                }
                if((i + N - 1) % 2 == color) {
                    //Right edge
                    ind = N * (i+1) - 1;
                    u[ind] = (h_sqr * f[ind] + u[ind-1] + u[ind+N] + u[ind-N])\
                            / 4;
                }
                for(j = 2 - (i+color)%2; j < N-1; j+=2) {
                    ind = N*i + j;
                    u[ind] = (h_sqr * f[ind] + u[ind-1] + u[ind+1] + u[ind+N]\
                                + u[ind-N]) / 4;
                }
            }
        }
    }
}
