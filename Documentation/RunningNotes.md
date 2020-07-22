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

#### 07-11-2019 ####
DM and JW corrected opt_alpha() as it was performing the "i in range(M)" loop in the particle dimension instead of the sample dimension

#### 14-11-2019 ####
# Start of document

#### 14-11-2019 ####
DM changed the max_L() function calls in opt_alpha() to take the phi0 output as phi0_guess argument

#### 6-5-2020 #####
Delayed update.

1) DM has changed the calc_K() function to use a Numpy pre-allocated matrix for the Kvals when memory allows it. Otherwise uses a block method which Numpy pre-allocates as far as possible before converting to sparse and then repeating the process, in the end stacking sparse matrices.

2) Commented a line in max_L() when "if phi0_guess == []" was True that added random noise to the phi guess.

3) Currently opt_alpha() uses DB89 sigma & mu when calling max_L() but deprojection_exe.py uses an initial guess for sigma & mu. Added option to specify whether deprojection_exe.py should use vars.ini guess values or not in vars.ini+deprojection_exe.py.

4) JW previously created a noniso version of calc_sigma2(), added option to specify noniso=1/0 (True or False) to both vars.ini and alpha_vars.ini, also became an argument for opt_alpha() and is pre-defined variable set to False in any calls to max_L()

5) sanity_check() printed output clarified and commented unnecessary lines. The function also claimed to be for comparing a pseudosamples gaussian distribution from DB89 with an estimated one. Changed description to specify that we compare real data DB89 values with deprojected ones.

6) Created a logfile that opt_alpha.py writes to, saving the final value and used values in alpha_vars.ini

#### 7-5-2020 #####
1) Added number of opt_alpha iteratios, fmin_cg calls and fmin_cg iterations to printed output of opt_alpha(). To do this, max_L() was modified to change its output based on which function/script is calling it.

2) Made a RUNS folder to which the Deproject_exe.py scipt automatically creates a folder for a run and saves the mxl data, and the plots. Also added a logfile that Deproject_exe.py writes to, saving the folder name for a run and the values which were used in vars.ini.

3) In the above logfile, also added the printed values of sanity_check()

#### 13-05-2020 ####

#### 14-05-2020 ####
- After discussing with John, DM determined that disp_guess + v_guess is only useful if we have an idea of what they should be, as this will create a phi0 that is closer to the real value and minimization will converge quicker.

#### 22-07-2020 ####
1) DM Added box_extrapolate to replace edge zeroes with extrapolated values in sec_der(). Noticed the same operation occurs in grad_sec_der() so implemented the extrapolation there as well.

2) DM made it so phi_guess now uses a combination of two Gaussians. This is because a single Gaussian cannot effectively be used with sec_der().

3) DM changed Kvals_selector so it will now use 90% of available RAM to operate.

4) DM wants to implement a multigrid approach, where phi_guess is the mxl calculated for a coarse grid and the grid is iteratively upsized. This would improve performance.
