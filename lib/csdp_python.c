#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "declarations.h"

double solve_sdp_python(k, block_num, block_sizes, a, rows, mat_inds, mat_vals)
int k;
int block_num;
int *block_sizes;
double *a;
int rows;
int *mat_inds;
double *mat_vals;
{
	int ret;
	int n;
	struct blockmatrix C;
	struct constraintmatrix *constraints;
	struct blockmatrix X,Z;
	double *y;
	double pobj,dobj;
	
	
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

	free_prob(n,k,C,a,constraints,X,y,Z);

	return(pobj);
}