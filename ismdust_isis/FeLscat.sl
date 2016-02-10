#! /usr/bin/env isis

%% FeLscat.sl -- Reads in table model data for FeL edge and applies it to a spectrum.
%%    This is a user defined model for fitting with ISIS.
%%
%% 2015.10.06 - lia@space.mit.edu
%%---------------------------------------------------------------------

variable FILE = "FeL_ext_res.dat";

private define read_xsect( filename )
{
    variable wavel, xext;
    (wavel, xext) = readcol( filename, 2, 3);

    return struct{wavel=wavel, tau=xext};
}

static variable X_lo, X_hi, Value;
static variable tauscat = read_xsect(SFILE);

define FeLedge_fit (lo, hi, par, fun)
{
    variable nH   = par[0];
    variable norm = par[1];

    variable Angs = 0.5*(lo+hi);
    variable tau  = nH * fsil * interpol(Angs,tauscat.wavel,tauscat.tau);

    return (1.0 - exp(-tau)) * fun;
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

add_slang_function ("FeLedge", ["nH","norm"]);
set_function_category("FeLedge", ISIS_FUN_OPERATOR);
set_param_default_hook("FeLedge", "FeLedge_defaults");

%provide ("silscat");