This provides a way to call CSDP from python.

CSDP is a software package for solving semidefinite programming problems. CSDP uses BLAS and LAPACK subroutines. The algorithm is a predictor-corrector version of the primal-dual barrier method of Helmberg, Rendl, Vanderbei, and Wolkowicz.

The original CSDP files, including the license and the name of the original authors, can be found in the csdp folder.

Usage:
solve_raw(k, block_sizes, a, mat_inds, mat_vals)
k - number of constraints
block_sizes - list of blocks sizes
a - right hand side values of the SDP
mat_inds - flattened list of the C and A matrix values with the sparse format (mat number, block number, first coord, second coord)
mat_vals - list of the corresponding values