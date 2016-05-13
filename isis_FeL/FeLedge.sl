#! /usr/bin/env isis

%% FeLscat.sl -- Reads in table model data for FeL edge and applies it to a spectrum.
%%    This is a user defined model for fitting with ISIS.
%%
%% 2015.10.06 - lia@space.mit.edu
%%---------------------------------------------------------------------

variable fepath = getenv("FEPATH");
variable FILE = fepath + "/FeL_xsect.dat";

private define read_xsect( filename )
{
    variable wavel, xext;
    (wavel, xext) = readcol( filename, 1, 3);

    return struct{wavel=wavel, tau=xext};
}

static variable X_lo, X_hi, Value;
static variable tauext = read_xsect(FILE);

define FeLedge_fit (lo, hi, par, fun)
{
    variable md   = par[0];

    variable Angs = 0.5*(lo+hi);
    variable tau  = md * interpol(Angs,tauext.wavel,tauext.tau);

    return fun * exp(-tau);
}

define FeLedge_defaults(i)
{
    switch(i)
    { case 0:
    return ( 1.0, 0, 0, 1000 );
    }
    { case 1:
    return ( 1.0, 1, 0, 1000 );
    }
}

add_slang_function ("FeLedge", ["md"]);
set_function_category("FeLedge", ISIS_FUN_OPERATOR);
set_param_default_hook("FeLedge", "FeLedge_defaults");

%provide ("silscat");