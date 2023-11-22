import ctypes
import array

libsdp = ctypes.CDLL('../lib/libsdp.so')

ci = ctypes.c_int
cd = ctypes.c_double
cpi = ctypes.POINTER(ci)
cpd = ctypes.POINTER(cd)


solver = libsdp.solve_sdp_python
solver.argtypes = [ci, ci, cpi, cpd, ci, cpi, cpd]
solver.restype = cd

def solve_raw(k, block_sizes, a, mat_inds, mat_vals):
	c_k = ci(k)
	c_block_num = ci(len(block_sizes))
	c_rows = ci(len(mat_vals))
	
	c_block_sizes = ctypes.cast(array.array('i', block_sizes).buffer_info()[0], cpi)
	c_a = ctypes.cast(array.array('d', a).buffer_info()[0], cpd)
	c_mat_inds = ctypes.cast(array.array('i', mat_inds).buffer_info()[0], cpi)
	c_mat_vals = ctypes.cast(array.array('d', mat_vals).buffer_info()[0], cpd)
	
	return solver(c_k, c_block_num, c_block_sizes, c_a, c_rows, c_mat_inds, c_mat_vals)