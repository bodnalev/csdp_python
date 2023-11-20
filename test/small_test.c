#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "declarations.h"

int main(argc,argv)
     int argc;
     char *argv[];
{
	int ret;
	int n;
	struct blockmatrix C;
	struct constraintmatrix *constraints;
	struct blockmatrix X,Z;
	double *y;
	double pobj,dobj;
	
	int k = 3;
	int block_num = 2;
	int block_sizes[] = {2, -4};
	double a[] = {0.0, -2.0/3, -1.0/3, 0.0};
	int rows = 12;
	int mat_inds[] = {0, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 2, 3, 1, 2, 2, 1, 2, 2, 2, 1, 2, 1, 1, 2, 2, 3, 3, 2, 2, 1, 1, 3, 2, 4, 4, 3, 2, 1, 1};
	double mat_vals[] = {-1.0, 1.0/3, 1.0/3, 1.0/3, 1.0/3, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0};
	
	ret=from_sparse_data(k,block_num,block_sizes,rows,mat_inds,mat_vals, &n, &C, &constraints, 1);
	if (ret != 0)
	{
		printf("Something is wrong with the problem\n");
		exit(201);
	};
	
	initsoln(n,k,C,a,constraints,&X,&y,&Z);
	
	ret=easy_sdp(n,k,C,a,constraints,0.0,&X,&y,&Z,&pobj,&dobj);
	
	/*
	somehow communicate back X, y, pobj, dobj etc, perhaps with copying write_sol() 
	*/
	
	printf("Result is: %f\n", pobj);
	free_prob(n,k,C,a,constraints,X,y,Z);
	
	return(ret);
}
