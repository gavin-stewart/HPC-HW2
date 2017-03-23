CC=gcc
OMP_CC=${CC} -fopenmp -lm
SOLVED_C=${wildcard omp_solved*.c}
GRID_SOLVERS_C=jacobi2D-omp.c gs2D-omp.c
SOLVED=${filter-out ${wildcard omp_solved*} omp_solved*.c}
GRID_SOLVERS=jacobi2D-omp gs2D-omp

all: ${SOLVED} ${GRID_SOLVERS}

%.c: %
	${OMP_CC} %.c -o %

clean:
	rm -f ${SOLVED} ${GRID_SOLVERS}


