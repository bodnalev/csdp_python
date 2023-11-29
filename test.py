from csdpy import solve_sdp

def test_quick():
	block_sizes = [2, -4]
	a = [-2.0/3, -1.0/3, 0.0]
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
	mat_vals = [-1.0, 1.0/3, 1.0/3, 1.0/3, 1.0/3, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0]
	sol = solve_sdp(block_sizes, a, mat_inds, mat_vals)
	print(sol)

if __name__ == "__main__":
	test_quick()