20/6-18:

get_L for 10 000 GDR1 stars when calc_K has a for-loop:

	3.2 s 	for N = 8x8x8, dv = 50x50x50, vmin = (-200,-200,200)
	
	23.1 s	for N = 50x50x50, dv = 8x8x8, vmin = (-200,-200,-200)
	
get_L for 10 000 GDR1 stars when calc_K just use arrays:

	2.6 s	for N = 8x8x8, dv = 50x50x50, vmin = (-200,-200,-200)
	
	16.8 s	for N = 50x50x50, dv = 8x8x8, vmin = (-200,-200,-200)
	
120718 results:

	Best run so far

	N = 1500
	v0 = 0,0,0
	disp = 50,50,50
	n = 25x25x25
	dv = 16x16x16
	v0_guess = 2,1,-3
	disp_guess = 30,30,30
	
	Sample mean ===> [-0.25048585  0.65953818 -2.12676714]
	Sample velocity dispersion ===> [ 50.45938417  50.68601034  49.58213394]
	Computed mean ===> [-0.07797956  0.62970625 -2.71046899]
	Computed velocity dispersion ===> [ 48.46015048  47.12892852  46.47691713]
