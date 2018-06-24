#! /usr/bin/env python

"""
    make_xs_ext.py -- A generic python script for making the extinction
        cross-section for use in XSPEC model ISMdust.

    Modify the magic numbers at the beginning of this file to change the
        grain size distribution and other dust properties.

    To remake the silicate cross-sections:
        ./make_xs_ext sil
    To remake the graphite cross-sections:
        ./make_xs_ext gra

    To construct the master extinction file for XSPEC model, run each command
    above, individually, and then run:
        ./make_xs_ext combine
    This will make edge_files/xs_ext_grid.fits from the two existing files:
        edge_files/silicate_xs.fits and edge_files/graphite_xs.fits

    Created by: Lia Corrales (lia@space.mit.edu)
    2015.11.18
    2016.06.25 -- updated for use with eblur/newdust
"""

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from scipy.interpolate import interp1d

import datetime
import multiprocessing
import os
import sys

## Requires installation of github.com/eblur/newdust
import newdust

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

_myname = 'liac@umich.edu'
_outdir = 'edge_files/' # file output directory

##-------------- Do not change these values --------------------------##

_mdust  = 1.e-4 # g cm^-2
_na     = 50    # number of grain sizes to use in distribution

## Astrophysical constants (Carroll & Ostlie)
_c    = 2.99792458e10  # cm s^-1
_h    = 6.6260755e-27  # erg s
_keV  = 1.60217733e-9 # erg/keV
_angs = 1.e-8          # cm/angs
_hc   = (_h*_c) / (_keV*_angs) # keV angs


## Final energy grid to use
## NOTE: ismdust model will break if you change the number of elements in the grid
_FeK   = np.arange(1.5, 1.91, 0.005) # 5 eV resolution in Fe K region
_Atemp = np.arange(1.0,130.0, 0.005) # 5 mAngstrom resolution everywhere else
_Etemp = _hc / (_Atemp[::-1])   # keV

# Final energy grid (keV)
_EGRID = np.append(np.append(_Etemp[_Etemp < _FeK[0]], _FeK),
    _Etemp[_Etemp > _FeK[-1]])

_FINAL_FILE = 'xs_ext_grid.fits'

print("Creating a final energy grid of length %d" % (len(_EGRID)) )
# in the past: 25800
# now (2018.06.24): 28405

# Low resulotion grids for computing
# This is for computing the cross-sections, then later
# we will interpolate those cross-sections onto _EGRID
_egrid_lores = np.logspace(-1.3, 1.1, 100.0)
# Energy grids for particular edges
_ANGSTROMS_OK   = np.linspace(22.0, 28.0, 1200) # 5 mA resolution
_ANGSTROMS_FeL  = np.linspace(15.0, 21.0, 1200) # 5 mA resolution
_ANGSTROMS_MgSi = np.linspace(5.0, 11.0, 1200) # 5 mA resolution
_ANGSTROMS_FeK  = np.linspace(1.5, 1.9, 1200) # 3.3 mA resolution
_ANGSTROMS_CK   = np.linspace(35, 48, 2600) # 5 mA resolution

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

def _insert_xsect(lores_grid, edge_grid, lores_xsect, edge_xsect):
    emin = edge_grid[0]
    emax = edge_grid[-1]

    result = np.array([])
    result = np.append(result, lores_xsect[lores_grid < emin])
    result = np.append(result, edge_xsect)
    result = np.append(result, lores_xsect[lores_grid > emax])
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
    thdulist.writeto(filename, overwrite=clobber)
    return

def _tau_scat_E( E, params ):
    amin, amax, p, rho, mdust, gtype = params
    result = ss.Kappascat(E=E, dist=_dustspec(params), scatm=ss.makeScatmodel('Mie',gtype))
    return result.kappa[0] * mdust

def _tau_ext_E( E, params ):
    amin, amax, p, rho, mdust, gtype = parmas
    result = ss.Kappaext(E=E, dist=_dustspec(params), scatm=ss.makeScatmodel('Mie',gtype))
    return result.kappa[0] * mdust

##-------------- Compute silicate values --------------##

def silicate_xs():
    egrid_sil = np.copy(_egrid_lores)
    region_list = [_ANGSTROMS_OK, _ANGSTROMS_FeL, _ANGSTROMS_MgSi, _ANGSTROMS_FeK]
    for edge in region_list:
        egrid_sil = _insert_edge_grid(egrid_sil, _hc/edge[::-1])

    sil_params = [_amin_s, _amax_s, _p_s, _rho_s, _mdust, 'Silicate']
    print("Making Silicate cross section with\n\tamin=%.3f\n\tamax=%.3f\n\tp=%.2f\n\trho=%.2f" \
          % (_amin_s, _amax_s, _p_s, _rho_s) )
    print("Output will be sent to %s" % (_outdir+_silfile))

    sil_comp = newdust.graindist.composition.CmSilicate(rho=_rho_s)
    sil_gpop = newdust.SingleGrainPop('Powerlaw', sil_comp, 'Mie',
        md=_mdust, amin=_amin_s, amax=_amax_s, p=_p_s)

    # do the calculation piece-by-piece
    sil_sca_by_reg = []
    sil_ext_by_reg = []
    for reg in region_list:
        sil_gpop.calculate_ext(_hc/reg[::-1], unit='kev')
        sil_sca_by_reg.append(sil_gpop.tau_sca)
        sil_ext_by_reg.append(sil_gpop.tau_ext)

    sil_gpop.calculate_ext(_egrid_lores, unit='kev')
    # cobble together the final cross-sections
    for reg, xsect in zip(region_list, sil_sca_by_reg):
        sca_sil = _insert_xsect(_egrid_lores, _hc/reg[::-1], sil_gpop.tau_sca, xsect)

    for reg, xsect in zip(region_list, sil_ext_by_reg):
        ext_sil = _insert_xsect(_egrid_lores, _hc/reg[::-1], sil_gpop.tau_ext, xsect)

    _write_all_xs_fits(_outdir+_silfile, egrid_sil, ext_sil, sca_sil, sil_params)
    return

##-------------- Compute graphite values --------------##

def graphite_xs():
    egrid_gra = np.copy(_egrid_lores)
    egrid_CK  = _hc/_ANGSTROMS_CK[::-1]
    egrid_gra = _insert_edge_grid(egrid_gra, egrid_CK)

    gra_params = [_amin_g, _amax_g, _p_g, _rho_g, _mdust, 'Graphite']
    print("Making Graphite cross section with\n\tamin=%.3f\n\tamax=%.3f\n\tp=%.2f\n\trho=%.2f" \
          % (_amin_g, _amax_g, _p_g, _rho_g))
    print("Output will be sent to %s" % (_outdir+_grafile))

    # Have to do the calculation for each orientation
    gra_para = newdust.graindist.composition.CmGraphite(rho=_rho_g, orient='para')
    gra_perp = newdust.graindist.composition.CmGraphite(rho=_rho_g, orient='perp')

    gra_gpop_para = newdust.SingleGrainPop('Powerlaw', gra_para, 'Mie',
        md=_mdust, amin=_amin_g, amax=_amax_g, p=_p_g)
    gra_gpop_perp = newdust.SingleGrainPop('Powerlaw', gra_perp, 'Mie',
        md=_mdust, amin=_amin_g, amax=_amax_g, p=_p_g)

    # Follow the graphie grain assumption Draine
    #   1/3 parallel and 2/3 perpendicular
    orient_frac = [1./3., 2./3.]

    # do the calculation for each orientation
    gra_sca, gra_ext = 0.0, 0.0 # lores version
    gra_sca_CK, gra_ext_CK = 0.0, 0.0 # hires C K edge region
    for f, gp in zip(orient_frac, [gra_gpop_para, gra_gpop_perp]):
        gp.calculate_ext(_hc/_ANGSTROMS_CK[::-1], unit='kev')
        gra_sca_CK += f * gp.tau_sca
        gra_ext_CK += f * gp.tau_ext
        gp.calculate_ext(_egrid_lores, unit='kev')
        gra_sca += f * gp.tau_sca
        gra_ext += f * gp.tau_ext

    # cobble together the final cross-sections
    ext_gra = _insert_xsect(_egrid_lores, egrid_CK, gra_ext, gra_ext_CK)
    sca_gra = _insert_xsect(_egrid_lores, egrid_CK, gra_sca, gra_sca_CK)

    # smooth the xs behavior for energies > esmooth
    def _smooth_xs(esmooth, xs, pslope):
        ipow   = np.where(egrid_gra >= esmooth)[0] # closest value to the desired esmooth value
        result = np.copy(xs)
        result[ipow] = xs[ipow[0]] * np.power(egrid_gra[ipow]/egrid_gra[ipow[0]], pslope)
        return result

    if _smooth_graphite_xs:
        ESMOOTH, PSCA, PABS = 1.0, -2.0, -2.9 # determined by hand
        print("Smoothing Graphite cross section with\n\tp=%.2f (scattering)\n\tp=%.2f (absorption)" % (PSCA,PABS))
        new_sca_gra = _smooth_xs(ESMOOTH, sca_gra, PSCA)
        new_abs_gra = _smooth_xs(ESMOOTH, ext_gra-sca_gra, PABS)
        new_ext_gra = new_sca_gra + new_abs_gra
        _write_all_xs_fits(_outdir+_grafile, egrid_gra, new_ext_gra, new_sca_gra, gra_params)
    else:
        _write_all_xs_fits(_outdir+_grafile, egrid_gra, ext_gra, sca_gra, gra_params)

    return

##-------------- Combine both into one fits file --------------##

def make_xs_fits(clobber=True):
    sil = Xsect(_outdir+_silfile)
    gra = Xsect(_outdir+_grafile)

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
    thdulist.writeto(_outdir+_FINAL_FILE, overwrite=clobber)
    return

##-------------- Main file execution ---------------------------##

if __name__ == '__main__':
    args = sys.argv
    #print(os.environ['PYTHONPATH'])
    if 'sil' in args:
        #print('silicate on')
        silicate_xs()
    if 'gra' in args:
        #print('graphite on')
        graphite_xs()
    if 'combine' in args:
        #print('combine on')
        make_xs_fits()
