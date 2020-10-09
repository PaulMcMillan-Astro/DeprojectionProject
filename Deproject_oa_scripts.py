import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import astropy.units as u
import astropy.coordinates as coord
import scipy.stats as st
import os
import inspect
import string
import builtins
import time
import psutil
from termcolor import colored
from astropy.table import Table
from astropy.io import ascii
from scipy import sparse as scisp
from scipy.optimize import fmin_cg
from scipy.interpolate import interpn
from scipy.ndimage import zoom
from datetime import date
from decimal import Decimal
from alpha_debugger import *
from IPython.display import clear_output
from Deproject_v1_0 import *

#@profile
def multigrid_steps_oa(n):
    step = 5
    while any(np.round(n/(2**step)) < 10):
        step -= 1
    box_steps = (2**np.linspace(0,step,step+1).reshape(step+1,1) * np.round(n/(2**step))).astype(int)
    
    
    return box_steps.astype(int)

#@profile
def max_L_oa(alpha, pvals, rhatvals, vmin, dv, n, phi0_guess=[], v0_guess=[], disp_guess=[], noniso=False, printing=True):
    builtins.L     = []
    builtins.gradL = [] 
    
    dvx, dvy, dvz = dv
    nx, ny, nz = n
    N = len(pvals)
    
    sigma2, vmean = calc_sigma2(pvals, rhatvals, True, noniso=noniso)
    if np.size(v0_guess) == 0:
        v0_guess = vmean
    if np.size(disp_guess) == 0:
        sigma = np.sqrt(sigma2)
        disp_guess = sigma
    if np.size(phi0_guess) == 0:
        phi0 = phi_guess(v0_guess, disp_guess, vmin, dv,n)
    else:
        phi0 = phi0_guess
     
    Kvals_args = (pvals, rhatvals, vmin, dv, n, N, printing)
    Kvals = Kvals_function_selector_oa(Kvals_args)

    args = (Kvals, N, alpha, dv, n, sigma2)
    phi0r = np.ravel(phi0)  # fmin_cg only takes one-dimensional inputs for the initial guess
    
    builtins.L.append(-1*get_neg_L(phi0r,Kvals, N, alpha, dv, n, sigma2))
    builtins.gradL.append(np.linalg.norm(get_grad_neg_L(phi0r,Kvals, N, alpha, dv, n, sigma2)))
        
    mxl, fopt, fcalls, gcalls, flag, phi_all = fmin_cg(get_neg_L, phi0r, fprime=get_grad_neg_L, gtol=1e-6, args=args, retall=True, disp=False, full_output=True, callback=callback)

    fmin_it = np.shape(phi_all)[0] - 1
                
    builtins.n = n
    builtins.dv = dv 
    
    mxlnew = mxl.reshape(n)
    
    return mxlnew, fmin_it


#@profile
def multigrid_max_L_oa(alpha, pvals, rhatvals, vmin, dv, n, phi0_guess=[], v0_guess=[], disp_guess=[], noniso=False, printing=False):
    N = len(pvals)
    vmax = vmin + dv*n
    
    box_steps = multigrid_steps_oa(n)
    n = box_steps[0]
    dv = (vmax-vmin)/n
    
    sigma2, vmean = calc_sigma2(pvals, rhatvals, True, noniso=noniso)
    if np.size(v0_guess) == 0:
        v0_guess = vmean
    if np.size(disp_guess) == 0:
        sigma = np.sqrt(sigma2)
        disp_guess = sigma
    if np.size(phi0_guess) == 0:
        phi0 = phi_guess(v0_guess, disp_guess, vmin, dv,n)
    else:
        phi0 = phi0_guess
        
    builtins.L     = [[] for _ in range(len(box_steps))]
    builtins.gradL = [[] for _ in range(len(box_steps))]
    fmin_it = 0
    for grid_step,n in enumerate(box_steps): 
        builtins.grid_step = grid_step    
        dv = (vmax-vmin)/n
            
        Kvals_args = (pvals, rhatvals, vmin, dv, n, N, printing)
        Kvals = Kvals_function_selector_oa(Kvals_args)
        
        args = (Kvals, N, alpha, dv, n, sigma2)
        phi0r = np.ravel(phi0)
        
        builtins.L[grid_step].append(-1*get_neg_L(phi0r,Kvals, N, alpha, dv, n, sigma2))
        builtins.gradL[grid_step].append(np.linalg.norm(get_grad_neg_L(phi0r,Kvals, N, alpha, dv, n, sigma2)))

        mxl, fopt, fcalls, gcalls, flag, phi_all = fmin_cg(get_neg_L, phi0r, fprime=get_grad_neg_L, gtol=1e-6, args=args, retall=True, disp=False, full_output=True, callback=callback_mg)
        fmin_it += np.shape(phi_all)[0] - 1   

        if grid_step != len(box_steps)-1:
            phi0 = zoomed_mxl(mxl.reshape(n))

    builtins.n = n
    builtins.dv = dv

    mxlnew = mxl.reshape(n)
    
    return mxlnew, fmin_it


#@profile
def Kvals_function_selector_oa(args):
    pvals, rhatvals, vmin, dv, n, N, printing = args
    
    AvMem = psutil.virtual_memory().available*0.9
    
    MaxFloats = AvMem / 8  # Number of floats we can handle assuming 8 bytes per float
    Nblock = int(np.floor(MaxFloats / np.prod(n)))  # Largest possible block size
    MemReq = 8*N*np.prod(n)/1e9
    AvMem /= 1e9
    if printing == True:
        print('Allocated RAM: %.2f GB  \nRequired RAM : %.2f GB | Block size = %s' % (AvMem, MemReq, Nblock))
    print('Allocated RAM: %.2f GB  \nRequired RAM : %.2f GB | Block size = %s' % (AvMem, MemReq, Nblock))

    if Nblock > N:
            Kvals = KvalsNumpyMethod(pvals, rhatvals, vmin, dv, n)
    else:
            Kvals = KvalsBlockMethod(pvals, rhatvals, vmin, dv, n, Nblock)
    
    return Kvals


#@profile
def opt_alpha(stdscr, alpha0, M, N, sample, vmin, dv, n, opt_tol=0.01, mise_tol=0.01, noniso=False):
    """Function that finds the optimal value of alpha for a given sample of stars.
    Given an initial guess of alpha, alpha0, it will draw M samples of size N from the resulting
    distribution f(v) computed using max_L. For each sample we perform the maximisation scheme
    to find the f(v) distribution for range of 10 alpha values between some upper and lower bounds. 
    The mean integrated square error (MISE) is then computed for each alpha value. The value with the lowest 
    MISE is the alpha_opt value and a new range of 10 alpha values is taken with upper and lower bounds set by
    the alpha values above and below opt_alpha in the previous range. The iterator repeats until the difference 
    |alpha_opt_previous-alpha_opt| falls below the mise tolerance at which point the alpha_opt value is set as the new 
    alpha_0 initial guess. The process then starts from a new initial range centered on the new alpha0. The iteration 
    ends when the difference between |alpha_0_previous-alpha_opt| falls below the optimization tolerance
    
    alpha0: The initial guess of alpha
    M: The number of samples to compute the MISE for
    N: Number of stars in each sample
    pvals: Array with the tangential velocity vectors for all the  original stars in our sample
    rhatvals: Array with the unit vector for each sample star
    vmin: Vector indicating the anchor of v-space
    dv: The dims of each cell
    n: The dims of our box
    mise_tol: The desired logarithmic tolerance for the minimisation scheme. default value is 0.01.
    opt_tol: The desired logarithmic tolerance for the optimization scheme. default value is 0.01.
    
    """
    ### Preparing counters
    ti = time.time()

    opt_diff = 10
    opt_it = 0

    mise_diff = 10
    mise_it = 0
    mise_it_tot = 0

    fmin_calls = 0
    fmin_it = 0

    ###### Iterator below ######
    # Setting up the initial range and guess
    logalpha_0 = np.log10(alpha0)
    logalpha_opt = logalpha_0
    logalphavals = np.linspace(logalpha_0 - 1, logalpha_0 + 1, 10)  # The initial set of alpha values

    vv, rrind, cartcoords, coordinds, pvals, rhatvals = read_params(sample, vmin, dv,
                                                                    n)  # Generate the needed params from input

    while opt_diff >= opt_tol:
        if mise_diff < mise_tol or mise_it == 0:
            if 'phi0' in locals():
                phi0, _, its = max_L(alpha0, pvals, rhatvals, vmin, dv, n, phi0_guess=phi0, disp=0, noniso=noniso)
                fmin_calls += 1
                fmin_it += its
            else:
                # print("phi0 not defined, using phi_guess for first phi0")
                phi0, _, its = max_L(alpha0, pvals, rhatvals, vmin, dv, n, disp=0, noniso=noniso)
                fmin_calls += 1
                fmin_it += its
            stdscr.addstr(2, 1, '%s' % "".join(([' '] * 42)))
            stdscr.addstr(5, 1, '%s' % "".join(([' '] * 42)))
            stdscr.addstr(2, 1, str('opt it       : %s' % opt_it), curses.color_pair(0) | curses.A_BOLD)
            stdscr.addstr(5, 1, str('time elapsed : %s hrs' % (np.around((time.time() - ti) / 3600, 3))),
                          curses.color_pair(0) | curses.A_BOLD)
            stdscr.refresh()

        vxvals, vyvals, vzvals, smpcoords, fv0 = sample_pdf(phi0, M, N, rrind, cartcoords, coordinds, vv,
                                                            n)  # Generate M pseudosamples of size N from phi0

        ise = np.zeros((M, len(logalphavals)))  # container for the integrated square errors

        for i in range(M):

            pspvals, psrhatvals = pseudosample_pvals(vxvals, vyvals, vzvals, smpcoords, i)  # Get a pseudosample

            for j in range(len(logalphavals)):
                alpha = 10 ** (logalphavals[j])

                phi, phiall, its = max_L(alpha, pspvals, psrhatvals, vmin, dv, n, phi0_guess=phi0, disp=0,
                                         noniso=noniso)
                fmin_calls += 1
                fmin_it += its

                fv = np.exp(phi)

                ise[i][j] = np.sum((fv - fv0) ** 2)

        mise = np.mean(ise, axis=0)

        # Finds alpha_opt using the mise and sets up new range
        xrange_old = logalphavals
        yrange_old = mise
        logalpha_opt, logalphavals, mise_diff, step = tenstep_mise(logalphavals, mise, logalpha_opt)
        xrange = logalphavals
        ### For the debugger:
        t = (time.time() - ti) / (3600)
        mise_it_tot += 1
        mise_it += 1
        make_string(stdscr, logalpha_0, opt_it, mise_it, mise_it_tot, mise_diff, xrange, xrange_old, yrange_old, t,
                    step, 'tenstep')

        if mise_diff < mise_tol:
            opt_diff = abs(logalpha_opt - logalpha_0)

            logalpha_0 = logalpha_opt  # We set the optimised alpha value to be our new initial guess
            alpha0 = 10 ** (logalpha_0)
            logalphavals = np.linspace(logalpha_0 - 1, logalpha_0 + 1, 10)

            builtins.Nps = M

            opt_it += 1  # One opt iteration completed

    alpha_fin = 10 ** (logalpha_0)

    stdscr.addstr(0, 1, str('The optimal value for alpha is : %.2f' % alpha_fin), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(1, 1, str('To run took : %s opt alpha iterations' % opt_it), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(2, 1, str('            : %s tenstep iterations' % mise_it_tot), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(3, 1, str('            : %s fmin_cg calls' % fmin_calls), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(4, 1, str('            : %s fmin_cg iterations' % fmin_it), curses.color_pair(0) | curses.A_BOLD)
    stdscr.refresh()

    
    return alpha_fin


#@profile
def opt_alpha_ternary(stdscr, alpha0, M, N, sample, vmin, dv, n, opt_tol=0.01, mise_tol=0.01, noniso=False):
    """Function that finds the optimal value of alpha for a given sample of stars using a ternary search algorithm.
    Given an initial guess of alpha, alpha0, it will draw M samples of size N from the resulting
    distribution f(v) computed using max_L. For each sample we perform the maximisation scheme
    to find the f(v) distribution for the ternary search left, left_third, right_third, and right alpha
    values. The mean integrated square error (MISE) is then computed for each ternary alpha  value to select the
    new ternary range. 
    
    This process repeats until the difference between |left-right| falls below the mise tolerance at which point
    the center of the ternary bounds is the new initial guess alpha_0. The ternary search then starts over until
    it finds |alpha_0_previous - alpha_0| is below the optimization tolerance
    
    alpha0: The initial guess of alpha
    M: The number of samples to compute the MISE for
    N: Number of stars in each sample
    pvals: Array with the tangential velocity vectors for all the  original stars in our sample
    rhatvals: Array with the unit vector for each sample star
    vmin: Vector indicating the anchor of v-space
    dv: The dims of each cell
    n: The dims of our box
    mise_tol: The desired logarithmic tolerance for the minimisation scheme. default value is 0.01.
    opt_tol: The desired logarithmic tolerance for the optimization scheme. default value is 0.01.
    
    """
    tree_file = make_treefile()
    with open(tree_file, 'a') as txt:
        txt.write('# Array shape: opt_it x log(alpha) x MISE\n')
    ### Preparing counters
    ti = time.time()

    opt_diff = 10
    opt_it = 0

    mise_diff = 10
    mise_it = 0
    mise_it_tot = 0

    fmin_calls = 0
    fmin_it = 0

    builtins.Nps = M
    ###### Iterator below ######
    # Setting up the initial range and guess

    logalpha_0 = np.log10(alpha0)
    xrange = np.linspace(logalpha_0 - 1, logalpha_0 + 1, 4)
    yrange = np.array([.0, .0, .0, .0])
    ind = np.array([True, True, True, True])

    logalphavals = xrange
    builtins.min_in_range = False
    vv, rrind, cartcoords, coordinds, pvals, rhatvals = read_params(sample, vmin, dv,
                                                                    n)  # Generate the needed params from input

    while opt_diff >= opt_tol:
        if mise_diff < mise_tol or mise_it == 0:
            if 'phi0' in locals():
                phi0, _, its = max_L(alpha0, pvals, rhatvals, vmin, dv, n, phi0_guess=phi0, disp=0, noniso=noniso)
                fmin_calls += 1
                fmin_it += its
            else:
                # print("phi0 not defined, using phi_guess for first phi0")
                phi0, _, its = max_L(alpha0, pvals, rhatvals, vmin, dv, n, disp=0, noniso=noniso)
                fmin_calls += 1
                fmin_it += its
            stdscr.addstr(2, 1, '%s' % "".join(([' '] * 42)))
            stdscr.addstr(5, 1, '%s' % "".join(([' '] * 42)))
            stdscr.addstr(2, 1, str('opt it       : %s' % opt_it), curses.color_pair(0) | curses.A_BOLD)
            stdscr.addstr(5, 1, str('time elapsed : %s hrs' % (np.around((time.time() - ti) / 3600, 3))),
                          curses.color_pair(0) | curses.A_BOLD)
            stdscr.refresh()

        vxvals, vyvals, vzvals, smpcoords, fv0 = sample_pdf(phi0, M, N, rrind, cartcoords, coordinds, vv,
                                                            n)  # Generate M pseudosamples of size N from phi0

        ise = np.zeros((Nps, 4))  # container for the integrated square errors

        for i in range(Nps):

            pspvals, psrhatvals = pseudosample_pvals(vxvals, vyvals, vzvals, smpcoords, i)  # Get a pseudosample

            for j in np.where(ind == True)[0]:
                alpha = 10 ** (logalphavals[j])

                phi, phiall, its = max_L(alpha, pspvals, psrhatvals, vmin, dv, n, phi0_guess=phi0, disp=0,
                                         noniso=noniso)
                fmin_calls += 1
                fmin_it += its

                fv = np.exp(phi)
                ise[i][j] = np.sum((fv - fv0) ** 2)

        mise = np.mean(ise, axis=0)

        np.copyto(yrange, mise, where=ind)

        with open(tree_file, 'a') as txt:
            data = np.ones(4, 3)
            data[:, 0] = np.ones(4, 1) * opt_it
            data[:, 1] = xrange
            data[:, 2] = yrange
            np.savetxt(txt, data)

        xrange_old = xrange
        yrange_old = yrange
        logalphavals, mise_diff, xrange, yrange, ind, step = ternary_mise(xrange, yrange)  # Get the new ternary range

        ### For the debugger:
        t = (time.time() - ti) / (3600)
        mise_it_tot += 1
        mise_it += 1
        make_string(stdscr, logalpha_0, opt_it, mise_it, mise_it_tot, mise_diff, xrange, xrange_old, yrange_old, t,
                    step, 'ternary')

        if mise_diff < mise_tol:
            logalpha_opt = (xrange[0] + xrange[3]) / 2

            opt_diff = abs(logalpha_opt - logalpha_0)

            logalpha_0 = logalpha_opt
            alpha0 = 10 ** (logalpha_0)

            xrange = np.linspace(logalpha_0 - 1, logalpha_0 + 1, 4)
            yrange = np.array([.0, .0, .0, .0])
            ind = np.array([True, True, True, True])

            logalphavals = xrange
            builtins.min_in_range = False
            builtins.Nps = M
            opt_it += 1

    alpha_fin = 10 ** (logalpha_0)

    stdscr.addstr(0, 1, str('The optimal value for alpha is : %.2f' % alpha_fin), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(1, 1, str('To run took : %s opt alpha iterations' % opt_it), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(2, 1, str('            : %s ternary iterations' % mise_it_tot), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(3, 1, str('            : %s fmin_cg calls' % fmin_calls), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(4, 1, str('            : %s fmin_cg iterations' % fmin_it), curses.color_pair(0) | curses.A_BOLD)
    stdscr.refresh()

    
    return alpha_fin


#@profile
def opt_alpha_gss(stdscr, alpha0, M, N, sample, vmin, dv, n, opt_tol=0.01, mise_tol=0.01, noniso=False):
    """Function that finds the optimal value of alpha for a given sample of stars using a golden-section search
    algorithm. Given an initial guess of alpha, alpha0, it will draw M samples of size N from the resulting
    distribution f(v) computed using max_L. For each sample we perform the maximisation scheme
    to find the f(v) distribution for the golden-section four values. They are:
    
    x1: Leftmost corner
    x2: 2nd leftmost point
    x3: 2nd rightmost point
    x4: Rightmost point
    
    They are specified by the intervals x2-x1:x3-x2:x4-x3 having the intervals widths in ratio 2-g:2g-3:2-g, where g 
    is the golden ratio. Similarly to a ternary search. if f(x3) < f(x2) the minimum is within x2 to x4 with new range
    x1',x2',x3',x4' = x2, x3, ?, x4. If f(x3) > f(x2) the minimum is within x1 to x3 with new range 
    x1',x2',x3',x4' = x1, ?, x2, x3. We calculate a new x3' or x2' point and estimate f() for it. The new range is
    chosen such that the widths follow the previous ratio.
    
    The mean integrated square error (MISE) is our f() and the process repeats until the new range limit difference
    falls below the mise tolerance, i.e. either |x2 - x4| or |x1 - x3|. At that point the central value (x3 or x2 in
    the previous cases) is set as the new guess for alpha0. The process starts over until it finds that 
    |alpha_0_previous - alpha_0| is below the optimization tolerance when the search concludes.
    
    alpha0: The initial guess of alpha
    M: The number of samples to compute the MISE for
    N: Number of stars in each sample
    pvals: Array with the tangential velocity vectors for all the  original stars in our sample
    rhatvals: Array with the unit vector for each sample star
    vmin: Vector indicating the anchor of v-space
    dv: The dims of each cell
    n: The dims of our box
    mise_tol: The desired logarithmic tolerance for the minimisation scheme. default value is 0.01.
    opt_tol: The desired logarithmic tolerance for the optimization scheme. default value is 0.01.
    
    """
    tree_file = make_treefile()
    try:
        s = open(tree_file, 'r').read().split('\n')
    except FileNotFoundError:
        with open(tree_file, 'a') as txt:
            txt.write('# Array shape: opt_it x log(alpha) x MISE\n')
    else:
        os.remove(tree_file)
        with open(tree_file, 'a') as txt:
            txt.write('# Array shape: opt_it x log(alpha) x MISE\n')
        stdscr.addstr(0, 51, ("Found an incomplete tree_file, removed it"),
                      curses.color_pair(4) | curses.A_BOLD | curses.A_UNDERLINE)
        stdscr.refresh()
    ### Preparing counters
    ti = time.time()

    opt_diff = 10
    opt_it = 0

    mise_diff = 10
    mise_it = 0
    mise_it_tot = 0

    fmin_calls = 0
    fmin_it = 0

    builtins.Nps = M
    ###### Iterator below ######
    # Setting up the initial range and guess

    logalpha_0 = np.log10(alpha0)

    invphi = (np.sqrt(5) - 1) / 2  # 1 / phi
    invphi2 = (3 - np.sqrt(5)) / 2  # 1 / phi^2

    x1 = logalpha_0 - 1
    x4 = logalpha_0 + 1
    x2 = x1 + invphi2 * (x4 - x1)
    x3 = x1 + invphi * (x4 - x1)
    xrange = np.array([x1, x2, x3, x4])
    yrange = np.array([.0, .0, .0, .0])
    ind = np.array([True, True, True, True])

    logalphavals = xrange
    builtins.min_in_range = False
    vv, rrind, cartcoords, coordinds, pvals, rhatvals = read_params(sample, vmin, dv,
                                                                    n)  # Generate the needed params from input

    while opt_diff >= opt_tol:
        if mise_diff < mise_tol or mise_it == 0:
            if 'phi0' in locals():
                phi0, _, its = max_L(alpha0, pvals, rhatvals, vmin, dv, n, phi0_guess=phi0, disp=0, noniso=noniso)
                fmin_calls += 1
                fmin_it += its
            else:
                # print("phi0 not defined, using phi_guess for first phi0")
                phi0, _, its = max_L(alpha0, pvals, rhatvals, vmin, dv, n, disp=0, noniso=noniso)
                fmin_calls += 1
                fmin_it += its
            stdscr.addstr(2, 1, '%s' % "".join(([' '] * 42)))
            stdscr.addstr(5, 1, '%s' % "".join(([' '] * 42)))
            stdscr.addstr(2, 1, str('opt it       : %s' % opt_it), curses.color_pair(0) | curses.A_BOLD)
            stdscr.addstr(5, 1, str('time elapsed : %s hrs' % (np.around((time.time() - ti) / 3600, 3))),
                          curses.color_pair(0) | curses.A_BOLD)
            stdscr.refresh()
        vxvals, vyvals, vzvals, smpcoords, fv0 = sample_pdf(phi0, Nps, N, rrind, cartcoords, coordinds, vv,
                                                            n)  # Generate M pseudosamples of size N from phi0
        ise = np.zeros((Nps, 4))  # container for the integrated square errors
        for i in range(Nps):

            pspvals, psrhatvals = pseudosample_pvals(vxvals, vyvals, vzvals, smpcoords, i)  # Get a pseudosample

            for j in np.where(ind == True)[0]:
                alpha = 10 ** (logalphavals[j])

                phi, phiall, its = max_L(alpha, pspvals, psrhatvals, vmin, dv, n, phi0_guess=phi0, disp=0,
                                         noniso=noniso)
                fmin_calls += 1
                fmin_it += its

                fv = np.exp(phi)
                ise[i][j] = np.sum((fv - fv0) ** 2)

        mise = np.mean(ise, axis=0)

        try:
            np.copyto(yrange, mise, where=ind)
        except TypeError:
            raise Exception('\nyrange : ' + str(yrange) + ', ' + str(yrange.dtype) +
                            '\nmise   : ' + str(mise) + ', ' + str(mise.dtype) +
                            '\nind    : ' + str(ind) + ', ' + str(ind.dtype))

        with open(tree_file, 'a') as txt:
            data = np.ones((4, 3))
            data[:, 0] = np.ones(4) * opt_it
            data[:, 1] = xrange
            data[:, 2] = yrange
            np.savetxt(txt, data)

        xrange_old = xrange
        yrange_old = yrange
        logalphavals, mise_diff, xrange, yrange, ind, step = gss_mise(xrange, yrange)  # Get the new gss range

        ### For the debugger:
        t = (time.time() - ti) / (3600)
        mise_it_tot += 1
        mise_it += 1
        make_string(stdscr, logalpha_0, opt_it, mise_it, mise_it_tot, mise_diff, xrange, xrange_old, yrange_old, t,
                    step, 'gss')

        if mise_diff < mise_tol:
            logalpha_opt = (xrange[0] + xrange[3]) / 2

            opt_diff = abs(logalpha_opt - logalpha_0)

            logalpha_0 = logalpha_opt
            alpha0 = 10 ** (logalpha_0)

            x1 = logalpha_0 - 1
            x4 = logalpha_0 + 1
            x2 = x1 + invphi2 * (x4 - x1)
            x3 = x1 + invphi * (x4 - x1)
            xrange = np.array([x1, x2, x3, x4])
            yrange = np.array([.0, .0, .0, .0])
            ind = np.array([True, True, True, True])

            logalphavals = xrange
            builtins.min_in_range = False
            builtins.Nps = M
            opt_it += 1

    alpha_fin = 10 ** (logalpha_0)

    with open(tree_file, 'a') as txt:
        txt.write('#FINISHED\n')

    stdscr.erase()
    stdscr.addstr(0, 1, str('The optimal value for alpha is : %.2f' % alpha_fin), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(1, 1, str('To run took : %s opt alpha iterations' % opt_it), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(2, 1, str('            : %s gss iterations' % mise_it_tot), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(3, 1, str('            : %s fmin_cg calls' % fmin_calls), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(4, 1, str('            : %s fmin_cg iterations' % fmin_it), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(6, 1, "Press any key 3 times to log output and close program!",
                  curses.color_pair(2) | curses.A_BOLD | curses.A_BLINK)
    stdscr.addstr(7, 1, "Press any key 3 times to log output and close program!",
                  curses.color_pair(2) | curses.A_BOLD | curses.A_BLINK)
    stdscr.addstr(8, 1, "Press any key 3 times to log output and close program!",
                  curses.color_pair(2) | curses.A_BOLD | curses.A_BLINK)
    stdscr.refresh()
    stdscr.getkey()
    stdscr.getkey()
    stdscr.getkey()
    
    return alpha_fin


#@profile
def read_params(sample, vmin, dv, n):
    pvals, rhatvals = calc_p_rhat(sample)
    vxmin, vymin, vzmin = vmin
    dvx, dvy, dvz = dv
    nx, ny, nz = n
    # We want to compute the rhat values for every sample
    # So we just draw M sets of size N of the coordinates from the original
    # sample that we will use for each iteration of the sample.

    sample.set_representation_cls(coord.CartesianRepresentation)
    cartcoords = np.zeros((len(sample), 3))
    cartcoords[:, 0] = sample.u
    cartcoords[:, 1] = sample.v
    cartcoords[:, 2] = sample.w
    coordinds = np.linspace(0, len(cartcoords) - 1, len(cartcoords))

    vxmax, vymax, vzmax = vxmin + nx * dvx, vymin + ny * dvy, vzmin + nz * dvz

    vx_bins = np.arange(vxmin, vxmax + dvx, dvx)
    vy_bins = np.arange(vymin, vymax + dvy, dvy)  # Bin-edges
    vz_bins = np.arange(vzmin, vzmax + dvz, dvz)

    vxc = (vx_bins[1:] + vx_bins[:-1]) / 2
    vyc = (vy_bins[1:] + vy_bins[:-1]) / 2
    vzc = (vz_bins[1:] + vz_bins[:-1]) / 2

    vxx, vyy, vzz = np.meshgrid(vxc, vyc, vzc, indexing='ij')  # The centre values of each cell in our v-space box
    vv = (vxx, vyy, vzz)
    ind = np.indices((nx, ny, nz))

    rind = np.ravel_multi_index(ind, (nx, ny, nz))  # An array containing the 3D coordinates of our box

    rrind = np.ravel(rind)
    
    return vv, rrind, cartcoords, coordinds, pvals, rhatvals


#@profile
def sample_pdf(phi0, M, N, rrind, cartcoords, coordinds, vv, n):
    nx, ny, nz = n

    fv0 = np.exp(phi0)  # The pdf f(v) given our inital guess alpha0
    fv0s = np.sum(fv0)
    prob = np.ravel(fv0 / fv0s)  # The normalised probability
    np.random.seed(0)
    smp = np.random.choice(rrind, (M, N), p=prob)  # Creation of M samples of size n given f(v)
    smpcoordinds = np.random.choice(coordinds, (M, N)).astype(int)  # We also draw random positions for these stars

    smpcoords = cartcoords[smpcoordinds]

    smpx, smpy, smpz = np.asarray(np.unravel_index(smp, (nx, ny, nz)))

    vxx, vyy, vzz = vv
    vxvals = vxx[smpx, smpy, smpz]
    vyvals = vyy[smpx, smpy, smpz]
    vzvals = vzz[smpx, smpy, smpz]
    
    return vxvals, vyvals, vzvals, smpcoords, fv0


#@profile
def pseudosample_pvals(vxvals, vyvals, vzvals, smpcoords, i):
    smpvx, smpvy, smpvz = vxvals[i], vyvals[i], vzvals[i]  # For every pseudosample we get velocities

    coordx, coordy, coordz = smpcoords[i, :, 0], smpcoords[i, :, 1], smpcoords[i, :, 2]  # ... as well as coordinates

    psample = coord.Galactic(u=coordx * u.pc, v=coordy * u.pc, w=coordz * u.pc, U=smpvx * (u.km / u.s),
                             V=smpvy * (u.km / u.s), W=smpvz * (u.km / u.s),
                             representation_type=coord.CartesianRepresentation,
                             differential_type=coord.CartesianDifferential)

    psample.set_representation_cls(coord.SphericalRepresentation, coord.SphericalCosLatDifferential)

    pspvals, psrhatvals = calc_p_rhat(psample)
    
    return pspvals, psrhatvals


#@profile
def tenstep_mise(logalphavals, mise, logalpha_opt_former):
    '''Function that performs the ten step zoom iteration and returns the new alpha0 and alpha range in log.
    Since our initial guess is very broad and covers 10 orders of magnitude, we have to narrow it down.
    This is done by taking the alpha values to the right and left of our optimal value to define the new range.
    If the optimal value is at the edges of our array, then we take the missing upper or lower bound
    to be 2 magnitudes larger than alpha_opt
    
    We also return the mise_diff to see if we are satisfied with the minimum. If we are at the edge however,
    this is set to 10 simply to ensure the iteration continues.'''

    optind = np.argwhere(mise == np.amin(mise))[0][0]  # The index of the optimimal alpha
    logalpha_opt = logalphavals[optind]
    mise_diff = abs(logalpha_opt_former - logalpha_opt)  # Check the improvement in the minimum from former best guess

    if (mise[:-1] < mise[1:]).all():
        lower = logalpha_opt - 1
        upper = logalphavals[optind + 1]
        mise_diff = 10
        step = 'Down'
    elif (mise[:-1] > mise[1:]).all():
        lower = logalphavals[optind - 1]
        upper = logalpha_opt + 1
        mise_diff = 10
        step = 'Up'
    elif (mise[:optind] > mise[1:optind + 1]).all() and (mise[optind:-1] < mise[optind + 1:]).all():
        lower = logalphavals[optind - 1]
        upper = logalphavals[optind + 1]
        step = 'Zoom'
    else:
        M = Nps
        builtins.Nps = M + 5
        lower = logalphavals[-1]
        upper = logalphavals[0]
        step = 'Again'

    return logalpha_opt, np.linspace(lower, upper, 10), mise_diff, step


#@profile
def ternary_mise(logalphavals, mise):
    '''Function that performs the ternary iteration and returns the new alpha0 and alpha range in log. If the
    lower(upper) edge has the smallest function value, we step 1 in log down(up) from left(right). I.e 
    the new range will be, using the old range values, left-1 <-> left_third (right_third <-> right+1)'''

    class Tern:
        #@profile
        def __init__(self, mise, alpha):
            self.x = alpha
            self.y = mise

    left, left_third, right_third, right = [Tern(mise[i], logalphavals[i]) for i in range(4)]

    # This is if either edge has the smallest mise and we need to expand the range to ensure we have the minimum
    # in our range of left <-> right.
    if (left.y < left_third.y < right_third.y < right.y):
        if not (min_in_range):
            xrange = np.linspace(left.x - 1, left_third.x, 4)
            yrange = np.array([0, 0, 0, left_third.y])
            logalphavals = np.array(list(xrange[1:3]) + [0])
            step = 'Down'
            return logalphavals, 10, xrange, yrange, np.array([True, True, True, False]), step
    elif (left.y > left_third.y > right_third.y > right.y):
        if not (min_in_range):
            xrange = np.linspace(right_third.x, right.x + 1, 4)
            yrange = np.array([right_third.x, 0, 0, 0])
            logalphavals = np.array([0] + list(xrange[:1]))
            step = 'Up'
            return logalphavals, 10, xrange, yrange, np.array([False, True, True, True]), step

    # This is if we have the minmum in the range left <-> right. Script only enters here in this case.
    elif ((left_third.y < left.y) and (right_third.y < right.y)):
        builtins.min_in_range = True
        if left_third.y < right_third.y:
            xrange = np.linspace(left.x, right_third.x, 4)
            yrange = np.array([left.y, 0, 0, right_third.y])
            logalphavals = np.array([0] + list(xrange[1:3]) + [0])
            step = 'Left'
            return logalphavals, abs(left.x - right_third.x), xrange, yrange, np.array([False, True, True, False]), step

        else:
            xrange = np.linspace(left_third.x, right.x, 4)
            yrange = np.array([left_third.y, 0, 0, right.y])
            logalphavals = np.array([0] + list(xrange[1:3]) + [0])
            step = 'Right'
            return logalphavals, abs(left_third.x - right.x), xrange, yrange, np.array([False, True, True, False]), step

    # If we go here the y-values imply a non-unimodal function, to which we must increase our
    # number of pseudosamples M and recalculate the MISE of the current range.
    else:
        M = Nps
        builtins.Nps = M + 5
        xrange = np.array([left.x, left_third.x, right_third.x, right.x])
        yrange = np.array([0, 0, 0, 0])
        step = 'Again'
        logalphavals = xrange
        return logalphavals, abs(xrange[0] - xrange[3]), xrange, yrange, np.array([True, True, True, True]), step


#@profile
def gss_mise(xrange, yrange):
    '''Function that performs the Golden section iteration and returns which new alpha needs to be evaluated, as well
    as its location within the array. If x1(x4) has the smallest function value, we step 1 in log down(up) from 
    x1(x4). I.e the new range will be, using the old range values, x1-1 <-> x2 (x3 <-> x4+1). Calculates mise_diff.'''

    x1, x2, x3, x4 = xrange
    y1, y2, y3, y4 = yrange
    # This is if either edge has the smallest mise and we need to expand the range to ensure we have the minimum
    # in our range of x1 <-> x4.
    if (y1 < y2 < y3 < y4):
        if not (min_in_range):
            x1 = x1 - 1;
            x4 = x2
            xrange = np.array(
                [x1, (x1 + (3 - np.sqrt(5)) / 2 * (x4 - x1)), (x1 + (np.sqrt(5) - 1) / 2 * (x4 - x1)), x4])
            yrange = np.array([0, 0, 0, y2],dtype=float)
            logalphavals = np.array(list(xrange[1:]) + [0])
            step = 'Down'
            return logalphavals, 10, xrange, yrange, np.array([True, True, True, False]), step
    elif (y1 > y2 > y3 > y4):
        if not (min_in_range):
            x1 = x3;
            x4 = x4 + 1
            xrange = np.array(
                [x1, (x1 + (3 - np.sqrt(5)) / 2 * (x4 - x1)), (x1 + (np.sqrt(5) - 1) / 2 * (x4 - x1)), x4])
            yrange = np.array([y3, 0, 0, 0],dtype=float)
            logalphavals = np.array([0] + list(xrange[1:]))
            step = 'Up'
            return logalphavals, 10, xrange, yrange, np.array([False, True, True, True]), step

    # This is if we have the minmum in the range left <-> right. Script only enters here in this case.
    elif ((y2 < y1) and (y3 < y4)):
        builtins.min_in_range = True
        if y2 < y3:
            xrange = np.array([x1, (x1 + (3 - np.sqrt(5)) / 2 * (x4 - x1)), x2, x3])
            yrange = np.array([y1, 0, y2, y3],dtype=float)
            step = 'Left'
            logalphavals = np.array([0, xrange[1], 0, 0])
            return logalphavals, abs(xrange[0] - xrange[3]), xrange, yrange, np.array([False, True, False, False]), step

        else:
            xrange = np.array([x2, x3, (x1 + (np.sqrt(5) - 1) / 2 * (x4 - x1)), x4])
            yrange = np.array([y2, y3, 0, y4])
            step = 'Right'
            logalphavals = np.array([0, 0, xrange[2], 0],dtype=float)
            return logalphavals, abs(xrange[0] - xrange[3]), xrange, yrange, np.array([False, False, True, False]), step

    # If we go here the y-values imply a non-unimodal function, to which we must increase our
    # number of pseudosamples M and recalculate the MISE of the current range.
    else:
        M = Nps
        builtins.Nps = M + 5
        xrange = np.array([x1, x2, x3, x4])
        yrange = np.array([0, 0, 0, 0],dtype=float)
        step = 'Again'
        logalphavals = xrange
        return logalphavals, abs(xrange[0] - xrange[3]), xrange, yrange, np.array([True, True, True, True]), step
