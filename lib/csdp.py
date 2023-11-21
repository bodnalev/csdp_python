import ctypes
import numpy as np

libsdp = ctypes.CDLL('./libsdp.so')

solve_sdp_python = libsdp.solve_sdp_python
solve_sdp_python.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double)]
solve_sdp_python.restype = ctypes.c_double

def csdp(k, block_sizes, a, mat_inds, mat_vals):
	c_k = ctypes.c_int(k)
	c_block_num = ctypes.c_int(len(block_sizes))
	c_block_sizes = np.array(block_sizes, dtype=np.intc).ctypes.data_as(ctypes.POINTER(ctypes.c_int))
	c_a = np.array(a, dtype=np.double).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
	c_rows = ctypes.c_int(len(mat_vals))
	c_mat_inds = np.array(mat_inds, dtype=np.intc).ctypes.data_as(ctypes.POINTER(ctypes.c_int))
	c_mat_vals = np.array(mat_vals, dtype=np.double).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
	return solve_sdp_python(c_k, c_block_num, c_block_sizes, c_a, c_rows, c_mat_inds, c_mat_vals)