#! /usr/bin/env python

"""
    make_xs_ext.py -- A generic python script for making the extinction 
        cross-section for use in XSPEC model ISMdust.
    
    Modify the magic numbers at the beginning of this file to change the 
        grain size distribution and other dust properties.
    
    Created by: Lia Corrales (lia@space.mit.edu)
    2015.11.18
"""

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from scipy.integrate import interp1d

import datetime
import multiproccesing
import os

## Requires installation of github.com/eblur/dust
import dust
import sigma_scat as ss

##-------------- Magic numbers, feel free to change -------------------##

## Silicate grain properties
_amin_s = 0.005 
_amax_s = 0.3  # min and max grain radius [um]
_p_s    = 3.5  # power law slope for grain size distribution
_rho_s  = 3.8  # dust grain density [g cm^-3]
_silfile = 'silicate_xs.fits'

## Graphite grain properties
_amin_g = 0.005
_amax_g = 0.3  # min and max grain radius [um]
_p_g    = 3.5  # power law slope for grain size distribution
_rho_g  = 2.2  # dust grain density [g cm^-3]
_grafile = 'graphite_xs.fits'
_smooth_graphite_xs = True 
# ^ Change to False to remove power law approximation 
#   on high E side of graphite cross section

_myname = 'lia@space.mit.edu'
_outdir = 'edge_files/' # file output directory

##-------------- Do not change these values --------------------------##

_mdust  = 1.e-4 # g cm^-2
_na     = 50    # number of grain sizes to use in distribution

## Astrophysical constants (Carroll & Ostlie)
_c    = 2.99792458e10  # cm s^-1
_h    = 6.6260755e-27  # erg s
_keV  = 1.60217733e-15 # erg
_angs = 1.e-8          # cm
_hc   = (_h*_c) / (_keV*_angs) # keV angs

## Final energy grid to use 
## NOTE: ismdust model will break if you change the number of elements in the grid
_dangs = 0.005
_AGRID = np.arange(1.0,130.0,__dangs) # wavelength [angs]
_EGRID = _hc / (_AGRID[::-1])       # keV
_FINAL_FILE = 'xs_ext_grid.fits'

print("Creating a final energy grid of length %d" % (len(_EGRID)) )

_egrid_lores = np.logspace(0.05, 12.0, 100.0)

# Energy grids for particular edges
_ANGSTROMS_OK   = np.linspace(22.0, 28.0, 1200)
_ANGSTROMS_FeL  = np.linspace(15.0, 21.0, 1200)
_ANGSTROMS_MgSi = np.linspace(5.0, 11.0, 1200)
_ANGSTROMS_CK   = np.linspace(35, 48, 2600)

##-------------- Supporting structures and functions -----------------##

class Xsect(object):
    def __init__(self, filename):
        data = fits.open(filename)[1].data
        self.energy = data['energy']
        self.tauext = data['ext']
    
    def __call__(self, e):
        result = interp1d(self.energy, self.tauext)
        return result(e)

def _insert_edge_grid(lores_grid, edge_grid):
    emin = edge_grid[0]
    emax = edge_grid[-1]

    result = np.array([])
    result = np.append(result, lores_grid[lores_grid < emin])
    result = np.append(result, edge_grid)
    result = np.append(result, lores_grid[lores_grid > emax])
    return result

def _write_all_xs_fits(filename, egrid, xs_ext, xs_sca, params, clobber=True):
    
    amin, amax, p, rho, mdust, gtype = params
    
    col1 = fits.Column(name='energy', format='E', array=egrid)
    col2 = fits.Column(name='angstroms', format='E', array=_hc/egrid)
    col3 = fits.Column(name='ext', format='E', array=xs_ext)
    col4 = fits.Column(name='sca', format='E', array=xs_sca)
    col5 = fits.Column(name='abs', format='E', array=xs_ext-xs_sca)
    
    cols  = fits.ColDefs([col1,col2,col3,col4,col5])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    #tbhdu.writeto(filename)
    
    prihdr = fits.Header()
    prihdr['AMIN']  = "%.3f" % (amin)
    prihdr['AMAX']  = "%.3f" % (amax)
    prihdr['P']     = "%.3f" % (p)
    prihdr['RHO']   = "%.3f" % (rho)
    prihdr['MDUST'] = "%.2e" % (mdust)
    prihdr['GTYPE'] = gtype
    prihdr['COMMENT'] = "Created by %s on %s" % (_myname, datetime.date.today())
    prihdu = fits.PrimaryHDU(header=prihdr)
    
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(filename, clobber=clobber)
    return

def _dustspec(params):
    amin, amax, p, rho, mdust, gtype = params
    radii = dust.Dustdist(rad=np.linspace(amin,amax,_na), p=p, rho=rho)
    return dust.Dustspectrum(md=mdust, rad=radii)

def _tau_scat_E( E, params ):
    amin, amax, p, rho, mdust, gtype = params
    result = ss.Kappascat(E=E, dist=_dustspec(params), scatm=ss.makeScatmodel('Mie',gtype))
    return result.kappa[0] * mdust

def _tau_ext_E( E, params ):
    amin, amax, p, rho, mdust, gtype = parmas
    result = ss.Kappaext(E=E, dist=_dustspec(params), scatm=ss.makeScatmodel('Mie',gtype))
    return result.kappa[0] * mdust

##-------------- Compute silicate values --------------##

def silicate_xs( nproc=4 ):
    egrid_sil = np.copy(_egrid_lores)
    for edge in [_ANGSTROMS_OK, _ANGSTROMS_FeL, _ANGSTROMS_MgSi]:
        egrid_sil = _insert_edge_grid(egrid_sil, _hc/edge[::-1])
    
    sil_params = [_amin_s, _amax_s, _p_s, _rho_s, _mdust, 'Silicate']
    print("Making Silicate cross section with\n\tamin=%.3d\n\tamax=%.3d\n\tp=%.2d\n\trho=%.2d" % sil_params)
    print("Output will be sent to %s" % (_silfile))
    
    """
    def _tau_sca(E):
        result = ss.Kappascat(E=E, dist=_dustspec(sil_params), scatm=ss.makeScatmodel('Mie','Silicate'))
        return result.kappa[0] * _mdust
    
    def _tau_ext(E):
        result = ss.Kappaext(E=E, dist=_dustspec(sil_params), scatm=ss.makeScatmodel('Mie','Silicate'))
        return result.kappa[0] * _mdust
    
    pool = multiprocessing.Pool(processes=nproc)
    sca_sil = pool.map(_tau_sca, egrid_sil)
    ext_sil = pool.map(_tau_ext, egrid_sil)
    pool.close()"""
    
    Ksca_sil = ss.Kappascat(E=egrid_gra, dist=_dustspec(sil_params), scatm=ss.makeScatmodel('Mie','Silicate'))
    Kext_sil = ss.Kappaext(E=egrid_gra, dist=_dustspec(sil_params), scatm=ss.makeScatmodel('Mie','Silicate'))
    
    sca_sil = Ksca_sil.kappa * _mdust
    ext_sil = Kext_sil.kappa * _mdust
    
    _write_all_xs_fits(_outdir+_silfile, egrid_sil, ext_sil, sca_sil, sil_params)
    return

##-------------- Compute graphite values --------------##

def graphite_xs( nproc=4 ):
    egrid_gra = np.copy(_egrid_lores)
    egrid_gra = _insert_edge_grid(egrid_gra, _hc/_ANGSTROMS_OK[::-1])
    
    gra_params = [_amin_g, _amax_g, _p_g, _rho_g, _mdust, 'Graphite']
    print("Making Graphite cross section with\n\tamin=%.3d\n\tamax=%.3d\n\tp=%.2d\n\trho=%.2d" % gra_params)
    print("Output will be sent to %s" % (_grafile))
    
    """
    def _tau_sca(E):
        result = ss.Kappascat(E=E, dist=_dustspec(gra_params), scatm=ss.makeScatmodel('Mie','Graphite'))
        return result.kappa[0] * _mdust
    
    def _tau_ext(E):
        result = ss.Kappascat(E=E, dist=_dustspec(gra_params), scatm=ss.makeScatmodel('Mie','Graphite'))
        return result.kappa[0]] * _mdust
    
    pool = multiprocessing.Pool(processes=nproc)
    sca_gra = pool.map(_tau_sca, egrid_gra)
    ext_sil = pool.map(_tau_ext, egrid_gra)
    pool.close()"""
    
    Ksca_gra = ss.Kappascat(E=egrid_gra, dist=_dustspec(gra_params), scatm=ss.makeScatmodel('Mie','Graphite'))
    Kext_gra = ss.Kappaext(E=egrid_gra, dist=_dustspec(gra_params), scatm=ss.makeScatmodel('Mie','Graphite'))
    
    sca_gra = Ksca_gra.kappa * _mdust
    ext_gra = Kext_gra.kappa * _mdust
    
    # smooth the xs behavior for energies > esmooth
    def _smooth_xs(esmooth, xs, pslope):
        ipow   = np.where(egrid_gra >= esmooth)[0] # closest value to the desired esmooth value
        result = np.copy(xs)
        result[ipow] = xs[ipow[0]] * np.power(egrid_gra[ipow]/egrid_gra[ipow[0]])
        return result
    
    if _smooth_graphite_xs = True:
        ESMOOTH, PSCA, PABS = 1.0, -2.0, -2.9 # determined by hand
        print("Smoothing Graphite cross section with\n\tp=%.2d (scattering)\n\tp=%.2d (absorption)" % (PSCA,PABS))
        new_sca_gra = _smooth_xs(ESMOOTH, sca_gra, PSCA)
        new_abs_gra = _smooth_xs(ESMOOTH, ext_gra-sca_gra, PABS)
        new_ext_gra = new_sca_gra + new_abs_gra
        _write_all_xs_fits(_outdir+_grafile, egrid_gra, new_ext_gra, new_sca_gra, gra_params)
    else:
        _write_all_xs_fits(_outdir+_grafile, egrid_gra, ext_gra, sca_gra, gra_params)
    
    return

##-------------- Combine both into one fits file --------------##

def make_xs_fits(clobber=True):
    sil = Xscat(_outdir+_silfile)
    gra = Xscat(_outdir+_grafile)
    
    col1 = fits.Column(name='energy', format='E', array=_EGRID*1.e3) # units of eV
    col2 = fits.Column(name='sil_ext', format='E', array=sil(_EGRID))
    col3 = fits.Column(name='gra_ext', format='E', array=gra(_EGRID))
    
    cols  = fits.ColDefs([col1,col2,col3])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    #tbhdu.writeto(filename)
    
    prihdr = fits.Header()
    prihdr['SIL_FILE']  = "%s" % (_silfile)
    prihdr['GRA_FILE']  = "%s" % (_grafile)
    prihdr['COMMENT'] = "Created by %s on %s" % (_myname, datetime.date.today())
    prihdu = fits.PrimaryHDU(header=prihdr)
    
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(_outdir+_FINAL_FILE, clobber=clobber)
    return

##-------------- Main file execution ---------------------------##

if __name__ == '__main__':
    sil=False, gra=False, combine=False
    args = os.argv[1:]
    if 'sil' in args:
        print('silicate on')
        #silicate_xs()
    if 'gra' in args:
        print('graphite on')
        #graphite_xs()
    if 'combine' in args:
        print('combine on')
        #make_xs_fits()


