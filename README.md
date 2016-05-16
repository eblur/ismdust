# ismdust

Contains calculations and fitting codes (ISIS table models and XSPEC models) from Corrales et. al. (2016)

http://arxiv.org/abs/1602.01100

To download, use git to create an ismdust folder:

    git clone git@github.com:eblur/ismdust.git ismdust

## Setup for ISMdust

Enter the ismdust directory and start XSPEC

    XSPEC12> initpackage ismdust lmodel_ismdust.dat .
    XSPEC12> exit

Finally, set an environment variable to point to location of ismdust, e.g.

    export ISMDUSTROOT=/path/to/ismdust/


## Setup for Fe-L edge fits with ISIS (Interactive Spectral Interpretation System) models

Add a line to your .isisrc file

    add_to_isis_load_path("/path/to/ismdust");

Set an environment variable to point to the location of the Fe-L edge templates:

    export FEPATH=/path/to/ismdust/ismdust_isis/

When you want to invoke the model, use the require function in ISIS to load the model.

    isis> require("FeLedge");

To invoke the model:

    isis> fit_fun("FeLedge(1, powerlaw(1))");
