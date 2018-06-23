#! /usr/bin/env isis

%% olivineabs.sl -- Reads in table model data olivine absorption and applies it.
%%    This is a user defined model for fitting with ISIS.
%%    If using this model, please cite Corrales et al. (2016)
%%    http://adsabs.harvard.edu/abs/2016MNRAS.458.1345C
%%    and Rogantini et al. (2018)
%%    http://adsabs.harvard.edu/abs/2018A%26A...609A..22R
%%
%% 2018.06.23 - lia@astro.wisc.edu
%%---------------------------------------------------------------------

_traceback = 1;

variable ISMPATH = getenv("ISMDUSTROOT");
variable OLIV_FILE = ISMPATH + "/edge_files/olivine_abs.fits";

private define read_xsect( filename )
{
    variable wavel, xabs;
    (wavel, xabs) = fits_read_col(filename, "angstroms", "abs");

    return struct{wavel=wavel[[::-1]], tau_abs=xabs[[::-1]]};
}

static variable X_lo, X_hi, Value;
static variable olv_ext = read_xsect(OLIV_FILE);

define olivineabs_fit(lo, hi, par, fun)
{
    variable olv_md  = par[0];
    variable Angs = 0.5*(lo+hi);
    variable tau  = olv_md * interpol(Angs,olv_ext.wavel,olv_ext.tau_abs);

    return fun * exp(-tau);
}

define olivineabs_defaults(i)
{
    switch(i)
    { case 0: % olv_md
    return ( 1.0, 0, 0, 10000 );
    }
}

add_slang_function ("olivineabs", ["olv_md"]);
set_function_category("olivineabs", ISIS_FUN_OPERATOR);
set_param_default_hook("olivineabs", "olivineabs_defaults");
