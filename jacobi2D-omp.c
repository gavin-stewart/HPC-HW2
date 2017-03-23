

#ifdef _OPENMP
#include <omp.h>
#endif
#define _USE_MATH_DEFINES

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "2Dgrid-verification.c"
#include "util.h"

#define TOL 1e-4

#ifdef _OPENMP
#define RUN_TYPE "parallel"
#else
#define RUN_TYPE "serial"
#endif

void compute_jacobi_iterate_serial(double*, double*, double*, int);
void compute_jacobi_iterate_parallel(double*, double*, double*, int);
double* run_jacobi(int, int, double*);

int main(int argc, char **argv){
    if(argc != 3) {
        printf("Usage: %s gridsize max_iterations\n", argv[0]);
        return EXIT_SUCCESS;
    }
    int i,j;
    int N = atoi(argv[1]);
    int max_iter = atoi(argv[2]);
    double *f = malloc(N * N * sizeof(double));
    double elapsed;
    timestamp_type t_start, t_end;
    for(i = 0; i < N; i++) {
        for(j = 0; j < N; j++) {
            f[N*i+j] = 1;
        }
    }
    get_timestamp(&t_start);
    double* sol = run_jacobi(N, max_iter, f);
    get_timestamp(&t_end);
    elapsed = timestamp_diff_in_seconds(t_start, t_end);
    printf("Elapsed time was %.3f seconds on %dx%d gridpoints and  <=%d\
 iterations with a %s Jacobi method.\n", elapsed, N, N, max_iter, RUN_TYPE);
    free(sol);
    free(f);
    return EXIT_SUCCESS;
}

void compute_jacobi_iterate(double *u, double *u_new, double *f, int N) {
    #ifdef _OPENMP
    compute_jacobi_iterate_parallel(u, u_new, f, N);
    #else
    compute_jacobi_iterate_serial(u, u_new, f, N);
    #endif
}

void compute_jacobi_iterate_serial(double *u, double *u_new, double *f, int N) {
    int i,j, ind;
    double h = 1.0 / (N + 1);
    double h_sqr = h * h;
    u_new[0] = 0.25 * (h_sqr * f[0] + u[1] + u[N]);
    for(j = 1; j < N-1; j++) {
        u_new[j] = 0.25 * (h_sqr * f[j] + u[j-1] + u[j+1] + u[j+N]);
    }
    u_new[N-1] = 0.25 * (h_sqr * f[N-1] + u[N-2] + u[2*N-1]);

    for(i = 1; i < N - 1; i++) {
        u_new[N*i] = 0.25 * (h_sqr * f[N*i] + u[N * (i - 1)] + u[N * (i + 1)]\
                     + u[N * i + 1]);
        for(j = 1; j < N - 1; j++) {
            ind = N * i + j;
            u_new[ind] = 0.25 * (h_sqr * f[ind] + u[ind + 1] + u[ind - 1]\
                                + u[ind + N] + u[ind - N]);
        }
        u_new[N*(i + 1) - 1] = 0.25 * (h_sqr * f[N*(i + 1) - 1] \
                    + u[N*(i + 1) - 2]  + u[N*i - 1] + u[N*(i + 2) - 1]);
    }

    u_new[N*(N-1)] = 0.25 * (h_sqr * f[N*(N-1)] + u[N * (N-1) + 1] \
                            + u[N*(N-2)]);
    for(j = 1; j < N-1; j++) {
        ind = N*(N-1) + j;
        u_new[ind] = 0.25 * (h_sqr * f[ind] + u[ind-1] + u[ind+1] \
                            + u[ind - N]);
    }
    u_new[N*N-1] = 0.25 * (h_sqr * f[N*N-1] + u[N*N-2] + u[N*(N-1)-1]);
}

void compute_jacobi_iterate_parallel(double *u, double *u_new, double* f, int N) {
    int i,j,ind;
    double h = 1.0 / (N+1);
    double h_sqr = h * h;
    #pragma omp parallel private(i,j,ind) shared(u, u_new, h, h_sqr, N, f)
    {
        #pragma omp sections nowait
        {
            #pragma omp section
            {
                //First row
                u_new[0] = 0.25 * (h_sqr * f[0] + u[1] + u[N]);
                for(j = 1; j < N-1; j++) {
                    u_new[j] = 0.25 * (h_sqr * f[j] + u[j-1] + u[j+1] + u[j+N]);
                }
                u_new[N-1] = 0.25 * (h_sqr * f[N-1] + u[N-2] + u[2*N-1]);
            }

            #pragma omp section
            {
                //Last row
                ind = N * (N-1);
                u_new[ind] = 0.25 * (h_sqr * f[ind] + u[ind + 1] + u[ind - N]);
                for(j = 1; j < N-1; j++) {
                    ind = N*(N-1) + j;
                    u_new[ind] = 0.25 * (h_sqr * f[ind] + u[ind - 1] \
                                            + u[ind + 1] + u[ind - N]);
                }
                ind = N*N - 1;
                u_new[ind] = 0.25 * (h_sqr * f[ind] + u[ind - 1]\
                                + u[ind - N]);
            }
        }

        #pragma omp for
        for(i=1; i < N-1; i++) {
            ind = N*i;
            u_new[ind] = 0.25 * (h_sqr * f[ind] + u[ind - N] + u[ind + N] \
                                + u[ind + 1]);
            for(j=1; j<N-1; j++) {
                ind = N*i + j;
                u_new[ind] = 0.25 * (h_sqr * f[ind] + u[ind-1] \
                                    + u[ind+1] + u[ind-N] + u[ind+N]);
            }
            ind = N * (i + 1) - 1;
            u_new[ind] = 0.25 * (h_sqr * f[ind] + u[ind - 1] + u[ind - N]\
                                + u[ind + N]);
        }
    }
}

double* run_jacobi(int N, int max_iter, double *f) {
    int iter_num;
    double *u, *u_new;
    double rel_res, start_residual;

    /* Setup */
    u = calloc(N * N,  sizeof(double));
    u_new = malloc(N * N * sizeof(double));

    //Compute the initial residual
    start_residual = compute_fd_residual(N, u, f);
    rel_res = 1.0;
    //Iterate
    for(iter_num = 0; iter_num < max_iter && rel_res > 1e-4; iter_num++) {
        //u_new is computed in the process of getting the residual.
        compute_jacobi_iterate(u, u_new, f, N);
        double *swp_tmp = u;
        u = u_new;
        u_new = swp_tmp;
        rel_res = compute_fd_residual(N, u, f) / start_residual;
    }
    printf("Final relative residual: %f\n", rel_res);
    free(u_new);
    return u;
}
