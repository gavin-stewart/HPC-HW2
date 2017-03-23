CC=gcc
OMP_CC=gcc -fopenmp -lm
SOLVED_C=${wildcard omp_solved*.c}
GRID_SOLVERS_C=jacobi2D-omp.c gs2D-omp.c
SOLVED=${filter-out  ${SOLVED_C},${wildcard omp_solved*}}
GRID_SOLVERS=jacobi2D-omp gs2D-omp
EXEC=${SOLVED} ${GRID_SOLVERS}

all: ${SOLVED} ${GRID_SOLVERS}

${EXEC} : %: %.c
	${OMP_CC} %.c -o %

clean:
	rm -f ${SOLVED} ${GRID_SOLVERS}


