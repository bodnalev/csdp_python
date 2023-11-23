/*
  Build the problem from the arrays.  Return 0 if ok, 1 if 
  failure.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "declarations.h"

void countentry();
int addentry();

int from_sparse_data(k,nblocks,block_sizes,rows,mat_inds,mat_vals,pn,pC,pconstraints,printlevel)
int k; 
int nblocks; 
int *block_sizes; 
int rows;
int *mat_inds; 
double *mat_vals;

int *pn;
struct blockmatrix *pC; 
struct constraintmatrix **pconstraints;
int printlevel;
{
	struct constraintmatrix *myconstraints;
	int i,j,ii;
	int blksz;
	int blk;
	int matno;
	int blkno;
	int indexi;
	int indexj;
	double ent;
	int ret;
	struct sparseblock *p;
	int *isdiag;
	double *tempdiag;
	
	/*
	* Number of blocks.
	*/
	if (nblocks<=0)
	{ 
		if (printlevel >= 1)
			printf("Incorrect SDPA, number of blocks must be positive\n");
		return(1);
	}
	
	/*
	* Keep track of which blocks have off diagonal entries. 
	*/
	isdiag=(int *)malloc((nblocks+1)*sizeof(int));
	for (i=1; i<=nblocks; i++)
		isdiag[i]=1;
	
	/*
	* Allocate space for the C matrix.
	*/
	pC->nblocks=nblocks;
	pC->blocks=(struct blockrec *)malloc((nblocks+1)*sizeof(struct blockrec));
	if (pC->blocks == NULL)
	{
		if (printlevel >= 1)
			printf("Storage allocation failed!\n");
		exit(205);
	}
	
	/*
	* Allocate space for the constraints.
	*/
	myconstraints=(struct constraintmatrix *)malloc((k+1)*sizeof(struct constraintmatrix));
	if (myconstraints == NULL)
	{
		if (printlevel >= 1)
			printf("Storage allocation failed!\n");
		exit(205);
	};
	
	/*
	* Null out all pointers in constraints.
	*/
	for (i=1; i<=k; i++)
	{
		myconstraints[i].blocks=NULL;
	};
	
	/*
	* And read the block structure.
	*/
	*pn=0;
	for (blk=1; blk<=nblocks; blk++)
	{
		blksz=block_sizes[blk-1];
		*pn=*pn+abs(blksz);
		if (blksz < 0)
		{
			pC->blocks[blk].blocksize=abs(blksz);
			pC->blocks[blk].blockcategory=DIAG;
			pC->blocks[blk].data.vec=(double *)malloc((1+abs(blksz))*sizeof(double));
			if (pC->blocks[blk].data.vec == NULL)
			{
				if (printlevel >= 1)
					printf("Storage allocation failed!\n");
				exit(205);
			};
			for (i=1; i<=abs(blksz); i++)
				pC->blocks[blk].data.vec[i]=0.0;

		}
		else
		{
			if (blksz > 46340)
			{
				printf("This problem is too large to be solved in I32LP64 mode.\n");
				exit(206);
			};
			
			pC->blocks[blk].blocksize=abs(blksz);
			pC->blocks[blk].blockcategory=MATRIX;
			pC->blocks[blk].data.mat=(double *)malloc((blksz*blksz)*sizeof(double));
			if (pC->blocks[blk].data.mat == NULL)
			{
				if (printlevel >= 1)
					printf("Storage allocation failed!\n");
				exit(205);
			};
			for (j=1; j<=blksz; j++)
				for (i=1; i<=blksz; i++)
					pC->blocks[blk].data.mat[ijtok(i,j,blksz)]=0.0;
		};
	};
	
	/*
	*  Now, loop through the entries, 
	*  counting entries in the constraint matrices block by block.
	* mat_inds,mat_vals
	*/
	
	
	for (ii=0; ii<rows; ii++)
	{
		matno = mat_inds[4*ii];
		blkno = mat_inds[4*ii+1];
		indexi = mat_inds[4*ii+2];
		indexj = mat_inds[4*ii+3];
		ent = mat_vals[ii];
		
		/*
		* Check the validity of these values.
		*/
		
		if ((matno < 0) || (matno > k) ||
		(blkno<1) || (blkno>nblocks) ||
		(indexi < 1) || (indexi > pC->blocks[blkno].blocksize) ||
		(indexj < 1) || (indexj > pC->blocks[blkno].blocksize))
		{
			if (printlevel >= 1)
				printf("Incorect SDPA. Bad values in line: %d %d %d %d %e \n",matno,blkno,indexi,indexj,ent);
			free(isdiag);
			return(1);
		};
		
		if (matno != 0)
		{
			if (ent != 0.0)
				countentry(myconstraints,matno,blkno,pC->blocks[blkno].blocksize);
		}
		else
		{
		/*
		* An entry in C. ignore it for now.
		*/
		};
	}
	/*
	* Now, go through each of the blks in each of the constraint matrices,
	* and allocate space for the entries and indices.
	*/
	for (i=1; i<=k; i++)
	{
		p = myconstraints[i].blocks;
		while (p != NULL)
		{
			/*
			* allocate storage for the entries in this block of this constraint.
			*/
			p->entries=(double *)malloc((p->numentries+1)*sizeof(double));
			if (p->entries == NULL)
			{
				if (printlevel >= 1)
					printf("Storage allocation failed!\n");
				exit(205);
			};
			p->iindices=(int *)malloc((p->numentries+1)*sizeof(int));
			if (p->iindices == NULL)
			{
				if (printlevel >= 1)
					printf("Storage allocation failed!\n");
				exit(205);
			};
			p->jindices=(int *)malloc((p->numentries+1)*sizeof(int));
			if (p->jindices == NULL)
			{
				if (printlevel >= 1)
					printf("Storage allocation failed!\n");
				exit(205);
			};
			p->numentries=0;
			p=p->next;
		};
	};
	
	
	
	
	
	
	
	
	

	/*
	*  Fill in the actual data.
	*/

	zero_mat(*pC);

	/*
	* Now, read the actual entries.
	*/
	
	for (ii=0; ii<rows; ii++)
	{
		matno = mat_inds[4*ii];
		blkno = mat_inds[4*ii+1];
		indexi = mat_inds[4*ii+2];
		indexj = mat_inds[4*ii+3];
		ent = mat_vals[ii];
		
		/*
		* Mark this block as not diagonal if indexi!=indexj.
		*/
		if ((indexi != indexj)  && (ent != 0.0))
			isdiag[blkno]=0;

		if (matno != 0)
		{
			if (ent != 0.0)
			{
				ret=addentry(myconstraints,matno,blkno,indexi,indexj,ent);
				if (ret != 0)
				{
					if (printlevel >= 1)
					{
						printf("Incorrect SDPA problem. Duplicate entry.\n");
						printf("matno=%d\n",matno);
						printf("blkno=%d\n",blkno);
						printf("indexi=%d\n",indexi);
						printf("indexj=%d\n",indexj);
					};
					free(isdiag);
					return(1);
				};
			};
		}
		else
		{
			/*
			* An entry in C. 
			*/
			if (ent != 0.0)
			{
				blksz=pC->blocks[blkno].blocksize;
				if (pC->blocks[blkno].blockcategory == DIAG)
				{
					if (pC->blocks[blkno].data.vec[indexi] != 0.0)
					{
						/*
						* We've got a duplicate entry in C!
						*/
						if (printlevel >= 1)
						{
							printf("Incorrect SDPA problem. Duplicate entry.\n");
							printf("matno=%d\n",matno);
							printf("blkno=%d\n",blkno);
							printf("indexi=%d\n",indexi);
							printf("indexj=%d\n",indexj);
						};
						free(isdiag);
						return(1);
					}
					else
						pC->blocks[blkno].data.vec[indexi]=ent;
				}
				else
				{
					if (pC->blocks[blkno].data.mat[ijtok(indexi,indexj,blksz)] != 0.0)
					{
						/*
						* We've got a duplicate entry in C!
						*/
						if (printlevel >= 1)
						{
							printf("Incorrect SDPA problem. Duplicate entry.\n");
							printf("matno=%d\n",matno);
							printf("blkno=%d\n",blkno);
							printf("indexi=%d\n",indexi);
							printf("indexj=%d\n",indexj);
						};
						free(isdiag);
						return(1);                    
					}
					else
					{
						pC->blocks[blkno].data.mat[ijtok(indexi,indexj,blksz)]=ent;
						pC->blocks[blkno].data.mat[ijtok(indexj,indexi,blksz)]=ent;
					};
				};
			};
		};
	}

	for (i=1; i<=nblocks; i++)
	{
		if ((pC->blocks[i].blockcategory != DIAG) && 
		(isdiag[i]==1) && (pC->blocks[i].blocksize > 1))
		{
			/*
			* We have a hidden diagonal block!
			*/
			if (printlevel >= 2)
				printf("Block %d is actually diagonal.\n",i);
			blksz=pC->blocks[i].blocksize;
			tempdiag=(double *)malloc((blksz+1)*sizeof(double));
			for (j=1; j<=blksz; j++)
				tempdiag[j]=pC->blocks[i].data.mat[ijtok(j,j,blksz)];
			free(pC->blocks[i].data.mat);
			pC->blocks[i].data.vec=tempdiag;
			pC->blocks[i].blockcategory=DIAG;
		};
	};

	/*
	* If the printlevel is high, print out info on constraints and block
	* matrix structure.
	*/
	if (printlevel >= 3)
	{
		printf("Block matrix structure.\n");
		for (blk=1; blk<=pC->nblocks; blk++)
		{
			if (pC->blocks[blk].blockcategory == DIAG)
				printf("Block %d, DIAG, %d \n",blk,pC->blocks[blk].blocksize);
			if (pC->blocks[blk].blockcategory == MATRIX)
				printf("Block %d, MATRIX, %d \n",blk,pC->blocks[blk].blocksize);
		};
	};
	/*
	* Next, setup issparse and NULL out all nextbyblock pointers.
	*/
	for (i=1; i<=k; i++)
	{
		p=myconstraints[i].blocks;
		while (p != NULL)
		{
			p->nextbyblock=NULL;
			p=p->next;
		};
	};
	/*
	* Free unneeded memory.
	*/
	free(isdiag);
	/*
	*  Put back all the returned values.
	*/
	*pconstraints=myconstraints;
	return(0);
}



