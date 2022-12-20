import sparse
import dask.array as da
import pickle
import numpy as np
import psutil
from Deproject_v1_0 import get_neg_L, get_grad_neg_L, calc_K

def unpack_pkl(pklfile):
    
    file = open(pklfile, 'rb')
    args = pickle.load(file)
    file.close
    return args.values()

def get_Nblock(n):
    # We allow only 80% of the available RAM to be used to allow other processes to run simulataneously.
    AvMem = psutil.virtual_memory().available*0.8

    MaxFloats = AvMem / 8  # Number of floats we can handle assuming 8 bytes per float
    Nblock = int(np.floor(MaxFloats / np.prod(n)))  # Largest possible block size
    return Nblock

def kvals_block_method(pvals, rhat, vmin, dv, n, nblock):
    n_stars = nblock*3

    # Run kvals for for all the blocks and stacking them sparsely after n_stars.
    kvals = None
    for block in np.array_split(np.arange(n_stars), n_stars//nblock):
        kvals_block = np.zeros((block.size, n.prod()))
        for i in block:
            kvals_block[i] = np.ravel(calc_K(pvals[i], rhat[i], vmin, dv, n))

            kvals_coo = scisp.coo_matrix(kvals_block)
            kvals = scisp.vstack((kvals, kvals_coo))
            del kvals_block
            kvals_block = np.zeros(kvals_block_shape)

    kvals_block = np.zeros((last_block, np.prod(n)))
    for i in range(last_block):
        kvals_block[i] = np.ravel(calc_K(pvals[i], rhat[i], vmin, dv, n))

    # Performing the final stack
    kvals_coo = scisp.coo_matrix(kvals_block)
    kvals = scisp.vstack((kvals, kvals_coo))
    kvals = scisp.csc_matrix(kvals)

    return kvals



phi0, pvals, rhatvals, vmin, alpha, dv, n, N, sigma2, printing = unpack_pkl('args.pkl')

