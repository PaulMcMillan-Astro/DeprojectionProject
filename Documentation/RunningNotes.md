# Running notes for the summer (de)project

It seemed like a good idea to have some running notes here that we can all add to. Be sure to update to the latest version before writing, and upload your changes when you're done.

Basics:

We're looking at using the techniques pioneered by Dehnen in [Local stellar kinematics from HIPPARCOS data](https://ui.adsabs.harvard.edu/#abs/1998MNRAS.298..387D/abstract) and, particularly [The Distribution of Nearby Stars in Velocity Space Inferred from HIPPARCOS Data](https://ui.adsabs.harvard.edu/#abs/1998AJ....115.2384D/abstract) to study the kinematics of the Solar neighborhood as found in Gaia DR2.

Our first task was determining how to calculate K(k|l) (D98, eq. 28) efficiently.

The working approach is to determine the values v_r which correspond to crossing the planes v_x = v_x,min + i * h_x; v_y = v_y,min + i * h_y etc., then putting these into an ordered list. The differences between consecutive values are then an ordered list of K values, which can be assigned to grid points by asking what box an intermediate value of vr would provide.

In the estimation of the second derivative (D98, eq. 30), we opted for an array-based method where we sum the terms from the adjacent cells for each cell l simultaneously and then add the l-specific term later.  

15/6/18: Bugs hunted down - 1) rhat must be a unit vector. 2) Vertices at the edge of the big box can cause trouble, because if the line crosses them then potentially there is a segment of the line (of zero length), assigned to a bin outside the grid.

19/6/18: We solved a problem with the estimation of the second derivative, where we wanted to avoid using for loops going through each cell in our v-space. This was done by summing all relevant terms for each box l in three separate arrays, one for each dimension.

20/6/18: Focused on performance and enhanced the calc_K function by removing the for-loop it contained. Computed L_tilde for 10 000 GDR1 stars in the test. Trying to think of a way to calculate K for all stars simultaneously.

21/6/18: Max_L is now fully implemented and runs. However it only does so for 2 iterations and returns values that does not make a lot of sense. Will implement a velocity model next week to investigate further.

25/6/18: Found that only one value of K was passed on to the max_L function which led to some funky stuff. Looks as if the function runs properly now as it does not terminate after 2 iterations anymore. It performed around 290 iterations for a sample of 3500 GDR1 stars on a 8x8x8 grid.

26/6/18: The max_L function runs for around 250 iterations and then returns the warning "Desired error not necessarily achieved due to precision loss". From the scipy.optimize documentation it seems that it cannot converge as the 'gradient and/or function calls were not changing'. It could also emerge from a bug where we get negative values of sum(exp(phi)*K(k|l)) in the first term, which doesn't seem possible. Will investigate this further.

3/7/18: After fixing an issue where the assumed values for the phi_guess was not correctly computed, the maximization now runs properly and optimizes successfully. It can reproduce simple Gaussian distributions for a low amount of stars, e.g. N=100. Requires further testing however.

6/7/18: The maximization code runs very slowly due to the fmin_cg function which takes a lot of time to numerically estimate the gradient of L_tilde. We therefore look to write a function that computes said gradient.

10/7/18: Finally made the max_L function behave as intended by fixing issues with get_grad_negL function. Further, there was a sign error in the return of the phi-values that maximize L_tilde. Added a few python scripts that contain various functions used.

	Deproject_v0.py: Contains all the base functions necessary to run max_L
	Deproject_plots.py: Contains all functions that plot (how about that)
	Deproject_test.py: As you might have guessed, this script contains the test functions.
	
11/7/18: Improved the various python files by adding the ability to use an input file of the same format as vars.ini. Also added a script called alpha_tester.py which is used to plot L and check the validity of alpha values, i.e. if it maximizes for the right velocity distribution.

12/7/18: Managed to greatly improve speed of max_L by rewriting some functions that used np.stack. We know instead create arrays of correct sizes and then assign values.

Thoughts 22/6/18: (We should have put minutes/plan of action from the meeting on Wednesday in here too.) An option when we first start working on getting the real data down from the Gaia Archive (and finding ways to choose slices in the HR diagram) is to start by determining the velocity dispersions of these slices and seeing how they vary across the HR diagram (rather than jumping straight to the estimation of Phi for different slices).

(Idea for much, much later which occurs to PJM - we can derive probabilities that stars are part of a velocity substructure, even without the line-of-sight velocities! We have a 'prior' P(v) which we derived from the full distribution, and can use this to say what the probability is for a given stars.

Forget I said this for now - it's for way, way later.
)
