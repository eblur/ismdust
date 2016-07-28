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

static variable X_lo, X_hi, Value;
static variable sil_ext = read_xsect(SIL_FILE);
static variable gra_ext = read_xsect(GRA_FILE);

define ismdust_fit(lo, hi, par, fun)
{
    variable sil_md   = par[0];
    variable gra_md   = par[1];

    variable Angs = 0.5*(lo+hi);
    variable tau  = sil_md * interpol(Angs,sil_ext.wavel,sil_ext.tau) +
                    gra_md * interpol(Angs,gra_ext.wavel,gra_ext.tau);

    return fun * exp(-tau);
}

define ismdust_defaults(i)
{
    switch(i)
    { case 0: % sil_md
    return ( 0.6, 0, 0, 10000 );
    }
    { case 1: % gra_md
    return ( 0.4, 0, 0, 10000 );
    }
}

add_slang_function ("ismdust", ["sil_md","gra_md"]);
set_function_category("ismdust", ISIS_FUN_OPERATOR);
set_param_default_hook("ismdust", "ismdust_defaults");

%%------- Absorption only model

define ismdust_abs_fit(lo, hi, par, fun)
{
    variable sil_md   = par[0];
    variable gra_md   = par[1];

    variable Angs = 0.5*(lo+hi);
    variable tau  = sil_md * interpol(Angs,sil_ext.wavel,sil_ext.tau_abs) +
                    gra_md * interpol(Angs,gra_ext.wavel,gra_ext.tau_abs);

    return fun * exp(-tau);
}

define ismdust_abs_defaults(i)
{
    switch(i)
    { case 0: % sil_md
    return ( 0.6, 0, 0, 10000 );
    }
    { case 1: % gra_md
    return ( 0.4, 0, 0, 10000 );
    }
}

add_slang_function ("ismdust_abs", ["sil_md","gra_md"]);
set_function_category("ismdust_abs", ISIS_FUN_OPERATOR);
set_param_default_hook("ismdust_abs", "ismdust_abs_defaults");

%%------- Scattering only model

define ismdust_sca_fit(lo, hi, par, fun)
{
    variable sil_md   = par[0];
    variable gra_md   = par[1];

    variable Angs = 0.5*(lo+hi);
    variable tau  = sil_md * interpol(Angs,sil_ext.wavel,sil_ext.tau_sca) +
                    gra_md * interpol(Angs,gra_ext.wavel,gra_ext.tau_sca);

    return fun * exp(-tau);
}

define ismdust_sca_defaults(i)
{
    switch(i)
    { case 0: % sil_md
    return ( 0.6, 0, 0, 10000 );
    }
    { case 1: % gra_md
    return ( 0.4, 0, 0, 10000 );
    }
}

add_slang_function ("ismdust_sca", ["sil_md","gra_md"]);
set_function_category("ismdust_sca", ISIS_FUN_OPERATOR);
set_param_default_hook("ismdust_sca", "ismdust_abs_defaults");
