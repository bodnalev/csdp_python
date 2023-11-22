from csdpy import solve_raw

def test_quick():
	k = 3
	block_sizes = [2, -4]
	a = [-2/3, -1/3, 0]
	mat_inds = [0, 2, 1, 1, 
				1, 1, 1, 1, 
				1, 1, 2, 1, 
				2, 1, 2, 1, 
				2, 1, 2, 2, 
				3, 1, 2, 2, 
				1, 2, 2, 2, 
				1, 2, 1, 1, 
				2, 2, 3, 3, 
				2, 2, 1, 1, 
				3, 2, 4, 4, 
				3, 2, 1, 1]
	mat_vals = [-1, 1/3, 1/3, 1/3, 1/3, 1, 1, -1, 1, -1, 1, -1]
	sol = solve_raw(k, block_sizes, a, mat_inds, mat_vals)
	assert abs(sol+0.5) < 1e-7, "Small test works fine"