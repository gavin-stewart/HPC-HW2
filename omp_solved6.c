/******************************************************************************
* FILE: omp_bug6.c
* DESCRIPTION:
*   This program compiles and runs fine, but produces the wrong result.
*   Compare to omp_orphan.c.
* AUTHOR: Blaise Barney  6/05
* LAST REVISED: 06/30/05
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define VECLEN 100

float a[VECLEN], b[VECLEN];

inline float dotprod ()
{
int i, tid;
float inner_sum = 0;

tid = omp_get_thread_num();
#pragma omp parallel shared(inner_sum) private(tid, i)
{
tid = omp_get_thread_num();
#pragma omp for reduction(+:inner_sum)
  for (i=0; i < VECLEN; i++)
    {
    inner_sum = inner_sum + (a[i]*b[i]);
    printf("  tid= %d i=%d\n",tid,i);
    }
}
return inner_sum;
}


int main (int argc, char *argv[]) {
int i;
float sum;

for (i=0; i < VECLEN; i++)
  a[i] = b[i] = 1.0 * i;

//No need to start the parallel section before the function call.  
//Move all parallelism into the function.
sum = dotprod();

printf("Sum = %f\n",sum);

}

