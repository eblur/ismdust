#! /usr/bin/env isis

%% ismdust.sl -- Reads in table model data for dust extinction edges and apply it to a spectrum.
%%    This is a user defined model for fitting with ISIS.
%%    If using this model, please cite Corrales et al. (2016)
%%    http://adsabs.harvard.edu/abs/2016MNRAS.458.1345C
%%
%% 2016.07.28 - lia@space.mit.edu
%%---------------------------------------------------------------------

_traceback = 1;

variable ISMPATH  = getenv("ISMDUSTROOT");
variable SIL_FILE = ISMPATH + "/edge_files/silicate_xs.fits";
variable GRA_FILE = ISMPATH + "/edge_files/graphite_xs.fits";

private define read_xsect( filename )
{
    variable wavel, xext, xsca, xabs;
    (wavel, xext, xabs, xsca) = fits_read_col(filename, "angstroms", "ext", "abs", "sca");

    return struct{wavel=wavel[[::-1]], tau=xext[[::-1]], tau_abs=xabs[[::-1]], tau_sca=xsca[[::-1]]};
}

private variable X_lo, X_hi, Value;
private variable sil_ext = read_xsect(SIL_FILE);
private variable gra_ext = read_xsect(GRA_FILE);

private define param_default_structure( value, freeze, pmin, pmax,
   hmin, hmax, pstep, prstep)
{
   variable param_def = struct{value, freeze, min, max, hard_min,
                               hard_max, step, relstep};
   param_def.value=value;
   param_def.freeze=freeze;
   param_def.min=pmin;
   param_def.max=pmax;
   param_def.hard_min=hmin;
   param_def.hard_max=hmax;
   param_def.step=pstep;
   param_def.relstep=prstep;
   return param_def;
}

define ismdust_fit(lo, hi, par)
{
    variable sil_md   = par[0];
    variable gra_md   = par[1];

    variable Angs = 0.5*(lo+hi);
    variable tau  = sil_md * interpol(Angs,sil_ext.wavel,sil_ext.tau) +
                    gra_md * interpol(Angs,gra_ext.wavel,gra_ext.tau);

    return exp(-tau);
}

define ismdust_defaults(i)
{
    switch(i)
    { case 0: % sil_md
       return param_default_structure( 0.6, 0, 0, 100, 0, 10000, 
                                       1.e-3, 1.e-3 );
    }
    { case 1: % gra_md
       return param_default_structure( 0.4, 0, 0, 100, 0, 10000, 
                                       1.e-3, 1.e-3 );
    }
}

add_slang_function ("ismdust", ["sil_md [1e-4 g/cm^2]","gra_md [1e-4 g/cm^2]"]);
set_param_default_hook("ismdust", "ismdust_defaults");

%%------- Absorption only model

define ismdust_abs_fit(lo, hi, par)
{
    variable sil_md   = par[0];
    variable gra_md   = par[1];

    variable Angs = 0.5*(lo+hi);
    variable tau  = sil_md * interpol(Angs,sil_ext.wavel,sil_ext.tau_abs) +
                    gra_md * interpol(Angs,gra_ext.wavel,gra_ext.tau_abs);

    return exp(-tau);
}

add_slang_function ("ismdust_abs", ["sil_md [1e-4 g/cm^2]","gra_md [1e-4 g/cm^2]"]);
set_param_default_hook("ismdust_abs", "ismdust_defaults");

%%------- Scattering only model

define ismdust_sca_fit(lo, hi, par)
{
    variable sil_md   = par[0];
    variable gra_md   = par[1];

    variable Angs = 0.5*(lo+hi);
    variable tau  = sil_md * interpol(Angs,sil_ext.wavel,sil_ext.tau_sca) +
                    gra_md * interpol(Angs,gra_ext.wavel,gra_ext.tau_sca);

    return exp(-tau);
}

add_slang_function ("ismdust_sca", ["sil_md [1e-4 g/cm^2]","gra_md [1e-4 g/cm^2]"]);
set_param_default_hook("ismdust_sca", "ismdust_defaults");

%% ------ Scattering halo spectral model

define ismdust_halo_fit(lo, hi, par)
{
    variable sil_md   = par[0];
    variable gra_md   = par[1];

    variable Angs = 0.5*(lo+hi);
    variable tau_sca = sil_md * interpol(Angs,sil_ext.wavel,sil_ext.tau_sca) +
                       gra_md * interpol(Angs,gra_ext.wavel,gra_ext.tau_sca);
    variable tau_abs = sil_md * interpol(Angs,sil_ext.wavel,sil_ext.tau_abs) +
                       gra_md * interpol(Angs,gra_ext.wavel,gra_ext.tau_abs);

    return exp(-tau_abs) * (1.0 - exp(-tau_sca));
}

add_slang_function ("ismdust_halo", ["sil_md [1e-4 g/cm^2]","gra_md [1e-4 g/cm^2]"]);
set_param_default_hook("ismdust_halo", "ismdust_defaults");

