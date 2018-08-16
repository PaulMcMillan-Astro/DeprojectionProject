#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 09:46:35 2018

@author: David Hobbs, Lund Observatory
"""
import numpy as np
import math;

mas2deg = 1.0/(3600*1000)
deg2rad = np.pi/180.0
mas2rad = mas2deg * deg2rad
    
Aprime  = np.asarray([[-0.0548755604162154, -0.8734370902348850, -0.4838350155487132],
                      [+0.4941094278755837, -0.4448296299600112, +0.7469822444972189],
                      [-0.8676661490190047, -0.1980763734312015, +0.4559837761750669]])

# Units are in deg or mas/(year) converted to radians/(year)
def transformGalIcrs(         lAr,      bAr,      varpiAr,      mulStarAr,      mubAr,
                     deltalStarAr, deltabAr, deltaVarpiAr, deltaMulStarAr, deltaMubAr):
                     
    N = len(lAr);
    raAr             = np.zeros(N)   
    decAr            = np.zeros(N)
    parallaxAr       = np.zeros(N)   
    pmraAr           = np.zeros(N)
    pmdecAr          = np.zeros(N)
    ra_errorAr       = np.zeros(N)
    dec_errorAr      = np.zeros(N)
    parallax_errorAr = np.zeros(N)
    pmra_errorAr     = np.zeros(N)
    pmdec_errorAr    = np.zeros(N)
    
    for i in range(N):
        l            = lAr[i]*deg2rad         
        b            = bAr[i]*deg2rad     
        varpi        = varpiAr[i]*mas2rad     
        mulStar      = mulStarAr[i]*mas2rad  
        mub          = mubAr[i]*mas2rad
        deltalStar   = deltalStarAr[i]*mas2rad
        deltab       = deltabAr[i]*mas2rad
        deltaVarpi   = deltaVarpiAr[i]*mas2rad
        deltaMulStar = deltaMulStarAr[i]*mas2rad
        deltaMub     = deltaMubAr[i]*mas2rad
        
        
        # Construct rGal
        rGal = np.asarray([[np.cos(l)*np.cos(b)],
                           [np.sin(l)*np.cos(b)],
                                     [np.sin(b)]])
        
        # Construct rIcrs
        rIcrs = np.matmul(np.transpose(Aprime),rGal)
        
        
        # Get alpha and delta
        alpha = math.atan2(rIcrs[1], rIcrs[0])
        delta = math.atan2(rIcrs[2], np.sqrt(rIcrs[0]**2 + rIcrs[1]**2))
        
        # The transformation of the proper motion components
        pIcrs = np.asarray([[-np.sin(alpha)],
                            [ np.cos(alpha)],
                            [ 0.0]])
        qIcrs = np.asarray([[-np.cos(alpha)*np.sin(delta)],
                            [-np.sin(alpha)*np.sin(delta)],
                            [               np.cos(delta)]])
        
        pGal  = np.asarray([[-np.sin(l)],
                            [ np.cos(l)],
                            [ 0.0]])
        qGal  = np.asarray([[-np.cos(l)*np.sin(b)],
                            [-np.sin(l)*np.sin(b)],
                            [           np.cos(b)]])
    

        muGal =  pGal*mulStar + qGal*mub
        muIcrs  = np.matmul(np.transpose(Aprime),muGal)
        
        # Icrs proper motions
        muAlphaStar = np.matmul(np.transpose(pIcrs),muIcrs)
        muDelta     = np.matmul(np.transpose(qIcrs),muIcrs)
        
        
        # Error propagation       
        g = np.asarray([[deltalStar], 
                        [deltab], 
                        [deltaVarpi], 
                        [deltaMulStar], 
                        [deltaMub]]);
        
        gal  = np.column_stack((pGal, qGal))
        icrs = np.column_stack((pIcrs, qIcrs))
        G = np.matmul(np.transpose(gal), np.matmul(Aprime,icrs))
        
        # Jacobian
        J = np.asarray([[G[0][0], G[0][1], 0.0,     0.0,     0.0],
                        [G[1][0], G[1][1], 0.0,     0.0,     0.0],
                        [    0.0,     0.0, 1.0,     0.0,     0.0],
                        [    0.0,     0.0, 0.0, G[0][0], G[0][1]],
                        [    0.0,     0.0, 0.0, G[1][0], G[1][1]]])
             
        #Errors in ICRS
        e = np.matmul(np.transpose(J),g);
        
        deltaAlphaStar   = e[0]
        deltaDelta       = e[1]
        deltaVarpi       = e[2]
        deltaMuAlphaStar = e[3]
        deltaMuDelta     = e[4]
        
        raAr[i]             = alpha/deg2rad    
        decAr[i]            = delta/deg2rad
        parallaxAr[i]       = varpi/mas2rad
        pmraAr[i]           = muAlphaStar/mas2rad
        pmdecAr[i]          = muDelta/mas2rad
        ra_errorAr[i]       = deltaAlphaStar/mas2rad
        dec_errorAr[i]      = deltaDelta/mas2rad
        parallax_errorAr[i] = deltaVarpi/mas2rad
        pmra_errorAr[i]     = deltaMuAlphaStar/mas2rad
        pmdec_errorAr[i]    = deltaMuDelta/mas2rad
        
    return       raAr,       decAr,       parallaxAr,       pmraAr,       pmdecAr, \
           ra_errorAr, dec_errorAr, parallax_errorAr, pmra_errorAr, pmdec_errorAr;



# Units are in deg or mas/(year) converted to radians/(year)
def transformIcrsGal(      ra,       dec,       parallax,       pmra,       pmdec,
                     ra_error, dec_error, parallax_error, pmra_error, pmdec_error):
 
    
    N = len(ra);
    lAr            = np.zeros(N)   
    bAr            = np.zeros(N)
    varpiAr        = np.zeros(N)   
    mulStarAr      = np.zeros(N)
    mubAr          = np.zeros(N)
    deltalStarAr   = np.zeros(N)
    deltabAr       = np.zeros(N)
    deltaVarpiAr   = np.zeros(N)
    deltaMulStarAr = np.zeros(N)
    deltaMubAr     = np.zeros(N)
    for i in range(N):
        alpha            = ra[i]*deg2rad         
        delta            = dec[i]*deg2rad     
        varpi            = parallax[i]*mas2rad     
        muAlphaStar      = pmra[i]*mas2rad  
        muDelta          = pmdec[i]*mas2rad
        deltaAlphaStar   = ra_error[i]*mas2rad
        deltaDelta       = dec_error[i]*mas2rad
        deltaVarpi       = parallax_error[i]*mas2rad
        deltaMuAlphaStar = pmra_error[i]*mas2rad
        deltaMuDelta     = pmdec_error[i]*mas2rad
        
        # Construct rIcrs
        rIcrs = np.asarray([[np.cos(alpha)*np.cos(delta)],
                            [np.sin(alpha)*np.cos(delta)],
                                          [np.sin(delta)]])
        
        # Construct rGal
        rGal = np.matmul(Aprime,rIcrs)
        
        # Get l and b
        l = math.atan2(rGal[1], rGal[0])
        b = math.atan2(rGal[2], np.sqrt(rGal[0]**2 + rGal[1]**2))
        
        # The transformation of the proper motion components
        pIcrs = np.asarray([[-np.sin(alpha)],
                            [ np.cos(alpha)],
                            [ 0.0]])
        qIcrs = np.asarray([[-np.cos(alpha)*np.sin(delta)],
                            [-np.sin(alpha)*np.sin(delta)],
                            [               np.cos(delta)]])
        
        pGal  = np.asarray([[-np.sin(l)],
                            [ np.cos(l)],
                            [ 0.0]])
        qGal  = np.asarray([[-np.cos(l)*np.sin(b)],
                            [-np.sin(l)*np.sin(b)],
                            [           np.cos(b)]])
        
        muIcrs =  pIcrs*muAlphaStar + qIcrs*muDelta
        muGal  = np.matmul(Aprime,muIcrs)
        
        # Galactic proper motions
        mulStar = np.matmul(np.transpose(pGal),muGal)
        mub     = np.matmul(np.transpose(qGal),muGal)
        
        
        # Error propagation       
        e = np.asarray([[deltaAlphaStar], 
                        [deltaDelta], 
                        [deltaVarpi], 
                        [deltaMuAlphaStar], 
                        [deltaMuDelta]]);
        

        gal  = np.column_stack((pGal, qGal))
    
        icrs = np.column_stack((pIcrs, qIcrs))
        G = np.matmul(np.transpose(gal), np.matmul(Aprime,icrs))
        
        # Jacobian
        J = np.asarray([[G[0][0], G[0][1], 0.0,     0.0,     0.0],
                        [G[1][0], G[1][1], 0.0,     0.0,     0.0],
                        [    0.0,     0.0, 1.0,     0.0,     0.0],
                        [    0.0,     0.0, 0.0, G[0][0], G[0][1]],
                        [    0.0,     0.0, 0.0, G[1][0], G[1][1]]])
             
        # Errors in Galactic
        g = np.matmul(J,e);
        
        deltalStar   = g[0]
        deltab       = g[1]
        deltaVarpi   = g[2]
        deltaMulStar = g[3]
        deltaMub     = g[4]
        
        lAr[i]            = l/deg2rad    
        bAr[i]            = b/deg2rad
        varpiAr[i]        = varpi/mas2rad
        mulStarAr[i]      = mulStar/mas2rad
        mubAr[i]          = mub/mas2rad
        deltalStarAr[i]   = deltalStar/mas2rad
        deltabAr[i]       = deltab/mas2rad
        deltaVarpiAr[i]   = deltaVarpi/mas2rad
        deltaMulStarAr[i] = deltaMulStar/mas2rad
        deltaMubAr[i]     = deltaMub/mas2rad
    
    return          lAr,      bAr,      varpiAr,      mulStarAr,      mubAr, \
           deltalStarAr, deltabAr, deltaVarpiAr, deltaMulStarAr, deltaMubAr;



## Define a source           
#ra             = np.asarray([55.80726241840959, 55.8522518960627])
#ra_error       = np.asarray([-4.563628563448031, 25.278002094871102]) #np.asarray([0.0595820642595418, 0.16002741908297363])
#dec            = np.asarray([24.533330029571047, 24.511626702125785])
#dec_error      = np.asarray([-5.229529085942313, -24.514670765238925]) #np.asarray([0.0382404923331015, 0.1128635706330098])
#parallax       = np.asarray([2.9765684942812154, 2.5782582828927305])
#parallax_error = np.asarray([2.9765684942812154, 2.5782582828927305]) #np.asarray([0.06800070746125245, 0.17666579992490822])
#pmra           = np.asarray([-4.563628563448031, 25.278002094871102])
#pmra_error     = np.asarray([-4.563628563448031, 25.278002094871102]) #np.asarray([0.11310783491759212, 0.33271813596762473])
#pmdec          = np.asarray([-5.229529085942313, -24.514670765238925])
#pmdec_error    = np.asarray([-5.229529085942313, -24.514670765238925]) #np.asarray([0.08850201342183286, 0.25382333995668005])
#
##ra             = np.asarray([000.00091185, 000.00379737])
##dec            = np.asarray([+01.08901332, -19.49883745])
##parallax       = np.asarray([3.54,          21.90])
##pmra           = np.asarray([-5.20,        181.21])
##pmdec          = np.asarray([-1.88,         -0.93])
##
##ra_error       = np.asarray([1.32,           1.28])
##dec_error      = np.asarray([0.74,           0.70])
##parallax_error = np.asarray([1.39,           3.10])
##pmra_error     = np.asarray([1.36,           1.74])
##pmdec_error    = np.asarray([0.81,           0.92])
#
#
##1	1	9.10	000.00091185	+01.08901332	3.54      -5.20 -1.88 1.32 0.74 1.39 1.36 0.81
##2	2	9.27	000.00379737	-19.49883745	21.90    181.21 -0.93 1.28 0.70 3.10 1.74 0.92
#
#
#print(ra, dec, parallax, pmra, pmdec)
#print(ra_error, dec_error, parallax_error, pmra_error, pmdec_error)
#
## Convert to Galactic
#l, b, varpi, mulStar, mub, deltalStar, deltab, deltaVarpi, deltaMulStar, deltaMub = \
#             transformIcrsGal(ra, dec, parallax, pmra, pmdec, ra_error, dec_error, parallax_error, pmra_error, pmdec_error);
#                 
#
#print(l, b, parallax, mulStar, mub)
#print(deltalStar, deltab, deltaVarpi, deltaMulStar, deltaMub)
#
## Check conversion with GalPy
#hc=coordinate.ICRS(ra=ra*u.degree,
#                   dec=dec*u.degree,
#                   pm_ra_cosdec=pmra*u.mas/u.yr,
#                   pm_dec=pmdec*u.mas/u.yr,)
#hgc=hc.transform_to(coordinate.Galactic)
#
#print(np.asarray(hgc.l), np.asarray(hgc.b), parallax, np.asarray(hgc.pm_l_cosb), np.asarray(hgc.pm_b))
#
## Convert back to ICRS
#ra, dec, parallax, pmra, pmdec, ra_error, dec_error, parallax_error, pmra_error, pmdec_error = \
#             transformGalIcrs(l, b, varpi, mulStar, mub, deltalStar, deltab, deltaVarpi, deltaMulStar, deltaMub);
#
#
#print(ra, dec, parallax, pmra, pmdec)
#print(ra_error, dec_error, parallax_error, pmra_error, pmdec_error)


