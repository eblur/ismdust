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

private define read_xsect( filename )
{
    variable wavel, xabs;
    (wavel, xabs) = fits_read_col(filename, "angstroms", "abs");

    return struct{wavel=wavel[[::-1]], tau_abs=xabs[[::-1]]};
}

private variable X_lo, X_hi, Value;
private variable olv_ext = read_xsect(OLIV_FILE);

define olivineabs_fit(lo, hi, par)
{
    variable olv_md  = par[0];
    variable Angs = 0.5*(lo+hi);
    variable tau  = olv_md * interpol(Angs,olv_ext.wavel,olv_ext.tau_abs);

    return exp(-tau);
}

define olivineabs_defaults(i)
{
    switch(i)
    { case 0: % olv_md
       return param_default_structure( 1.0, 0, 0, 100, 0, 10000, 
                                       1.e-3, 1.e-3 );
    }
}

add_slang_function ("olivineabs", ["olv_md [1.e-4 g/cm^2]"]);
set_param_default_hook("olivineabs", "olivineabs_defaults");
