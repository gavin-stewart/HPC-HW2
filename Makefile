CC=gcc
CFLAGS=-fopenmp -lm
SOLVED_C=${wildcard omp_solved*.c}
GRID_SOLVERS_C=jacobi2D-omp.c gs2D-omp.c
SOLVED=${SOLVED_C:.c=}
GRID_SOLVERS=jacobi2D-omp gs2D-omp
EXEC=${SOLVED} ${GRID_SOLVERS}

all: ${EXEC}

${SOLVED} : % : %.c
	${CC} ${CFLAGS} $^ -o $@

${GRID_SOLVERS} : % : %.c
ifeq (${SERIAL}, [tT])
	${CC} ${CFLAGS} $^ -o $@
else
	${CC} -lm $^ -o $@
endif

clean:
	rm -f ${SOLVED} ${GRID_SOLVERS}


