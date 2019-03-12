#! /usr/bin/env isis

%% FeLscat.sl -- Reads in table model data for FeL edge and applies it to a spectrum.
%%    This is a user defined model for fitting with ISIS.
%%
%% 2015.10.06 - lia@space.mit.edu
%%---------------------------------------------------------------------

variable fepath = getenv("FEPATH");
variable FILE = fepath + "/FeL_xsect.dat";

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
    (wavel, xabs) = readcol( filename, 1, 5);

    return struct{wavel=wavel, tau=xabs};
}

private variable X_lo, X_hi, Value;
private variable tauabs = read_xsect(FILE);

define FeLabs_fit (lo, hi, par)
{
    variable md   = par[0];

    variable Angs = 0.5*(lo+hi);
    variable tau  = md * interpol(Angs,tauabs.wavel,tauabs.tau);

    return exp(-tau);
}

define FeLabs_defaults(i)
{
    switch(i)
    { case 0:
       return param_default_structure( 1.0, 0, 0, 100, 0, 10000, 
                                       1.e-3, 1.e-3 );
    }
    { case 1:
       return param_default_structure( 1.0, 1, 0, 100, 0, 10000, 
                                       1.e-3, 1.e-3 );
    }
}

add_slang_function ("FeLabs", ["md [1e-4 g/cm^2]"]);
set_param_default_hook("FeLabs", "FeLabs_defaults");

