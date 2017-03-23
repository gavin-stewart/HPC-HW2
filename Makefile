CC=gcc
OMP_CC=gcc -fopenmp -lm
SOLVED_C=${wildcard omp_solved*.c}
GRID_SOLVERS_C=jacobi2D-omp.c gs2D-omp.c
SOLVED=${SOLVED_C:.c=}
GRID_SOLVERS=jacobi2D-omp gs2D-omp
EXEC=${SOLVED} ${GRID_SOLVERS}

all: ${EXEC}

${EXEC} : % : %.c
	${OMP_CC} $^ -o $@

clean:
	rm -f ${SOLVED} ${GRID_SOLVERS}


