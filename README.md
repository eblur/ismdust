# ismdust

Contains calculations and fitting codes (ISIS table models and XSPEC models) from Corrales et. al. (2016)

http://arxiv.org/abs/1602.01100

To download, use git to create an ismdust folder:

        git clone https://github.com/eblur/ismdust.git ismdust

## About the model parameters

Each model applies the value `exp(-tau md)` to the flux, where tau is the path-integrated
absorption, scattering, or extinction optical depth and `md` is the
dust mass column in units of 10^-4 g cm^-2.

To estimate the appropriate dust mass column for a given ISM column (NH),
one must choose a dust-to-gas mass ratio. For for the Milky Way ISM,
this value is typically ~0.01.

For example, NH = 10^22 cm^-2 corresponds to NH * m_p (proton mass) * 0.01 =
 1.67e-4 g cm^-2. One can then set, for example, the `md_sil`
parameter in ISMdust to 1.67.

**ISMdust model parameters**

    sil_md : dust mass column for Silicate in units of 1.e-4
      default = 0.6

    gra_md : dust mass column for Graphite in units of 1.e-4
      default = 0.4

    redshift : redshift of the obscuring dust (XSPEC only)
      default = 0.0

The default mixture of 60%/40% silicate/graphite is described in
[Corrales et al. (2016)](http://arxiv.org/abs/1602.01100).
Please cite this paper if you use this model.

**olivineabs model parameters**

  olv_md : dust mass column for Olivine grains in units of 1.e-4
    default = 1.0

  redshift : redshift of the obscuring dust (XSPEC only)
    default = 0.0

A simple model for olivine absorption that uses the silicate absorption model
from ISMdust with the [Rogantini et al. (2018)](http://adsabs.harvard.edu/abs/2018A%26A...609A..22R)
cross-section for the Fe K edge incorporated.
Please cite both Corrales et al. (2016) and Rogantini et al. (2018)
if you use this model.

## XSPEC setup

Enter the ismdust directory and start XSPEC

    XSPEC12> initpackage ismdust lmodel_ismdust.dat .
    XSPEC12> exit

Finally, set an environment variable to point to location of ismdust, e.g.

    export ISMDUSTROOT=/path/to/ismdust/

Now, when you want to load the ISMDUST model in XSPEC, load the correct library. The following example loads ismdust and then applies the extinction model to a power law component:

    XSPEC12> lmod ismdust $::env(ISMDUSTROOT)
    XSPEC12> mo ismdust*pow

Note: You can omit the `$::env(ISMDUSTROOT)` portion if you have the `LMODDIR`
environment variable set to the location where ISMdust is installed.

Try the test file to make sure it's working.

    XSPEC12> @test.xcm

### For Silicate absorption with Olivine Fe K cross-section

Follow the instructions above to set your `ISMDUSTROOT` environment variable
(and `LMODDIR`, if you choose).

Now install the local model

      XSPEC12> initpackage olivineabs lmodel_olivineabs.dat .

To load the model in XSPEC:

    XSPEC12> lmod olivineabs $::env(ISMDUSTROOT)
    XSPEC12> mo olivineabs*pow

Try the test file to make sure it's working.

    XSPEC12> @test_olivine.xcm

## ISIS (Interactive Spectral Interpretation System) setup

Add a line to your .isisrc file

    add_to_isis_load_path("/path/to/ismdust/ismdust_isis");

Set an environment variable (same way you would for the XSPEC model)

    export ISMDUSTROOT=/path/to/ismdust/

When you want to invoke the model, use the require function in ISIS to load ismdust

    isis> require("ismdust");

To set up the model extinction model with a power law continuum, for example, do:

    isis> require("ismdust");
    isis> fit_fun("ismdust(1, powerlaw(1))");

**For absorption component only**

    isis> require("ismdust");
    isis> fit_fun("ismdust_abs(1, powerlaw(1))");

**For scattering component only**

    isis> require("ismdust");
    isis> fit_fun("ismdust_sca(1, powerlaw(1))");

See `ismdust_isis/test_ismdust.sl`

### For Silicate absorption with Olivine Fe K cross-section

Use the same set up instructions as above.

To run the model with a power law continuum, for example, do:

    isis> require("olivineabs");
    isis> fit_fun("olivineabs(1, powerlaw(1))");

See also `ismdust_isis/test_olivine.sl`

### Fe-L edge fits with ISIS

Set an environment variable to point to the location of the Fe-L edge templates:

    export FEPATH=/path/to/ismdust/ismdust_isis/

To invoke the model:

    isis> require("FeLedge");
    isis> fit_fun("FeLedge(1, powerlaw(1))");
