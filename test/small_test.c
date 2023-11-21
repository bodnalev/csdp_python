#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "declarations.h"

int main(argc,argv)
     int argc;
     char *argv[];
{
	int k = 3;
	int block_num = 2;
	int rows = 12;
	
	
	int block_sizes[] = {2, -4};
	double a[] = {-2.0/3, -1.0/3, 0.0};
	int mat_inds[] = {0, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 2, 3, 1, 2, 2, 1, 2, 2, 2, 1, 2, 1, 1, 2, 2, 3, 3, 2, 2, 1, 1, 3, 2, 4, 4, 3, 2, 1, 1};
	double mat_vals[] = {-1.0, 1.0/3, 1.0/3, 1.0/3, 1.0/3, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0};
	
	printf("Result is: %f\n", solve_sdp_python(k, block_num, block_sizes, a, rows, mat_inds, mat_vals));
	
	return(0);
}
