ICInG ( **I**solated Collapse Initial Conditions Generator)

---------------
Description
--------------

Script created by Michael Tremmel (University of Washington) to generate initial conditions for a cosmological collapsing, isolated halo. Can have any number of DM and Gas particles the user chooses.

Formation of the initial conditions follows Evrard 1988 and assumes an Einstein-deSitter cosmology.

Uses a glass distribution to avoid numerical artifacts

Includes the ability to impart angular momentum as well as random motions to help avoid radial instabilities

This was developed for use with Gasoline and ChaNGa N-body+SPH simulations, but should be able to be ported to other formats.

type help(pyCosmoColl.IC.MkCosmoSphereICs) in the python environment to get information on parameters to input.

-------------
A Note on Units
-------------

Uses the following units by default:

Mass Unit            = 1.84793e16 Msol

Distance Unit        = 50000 kpc.



rho_crit is assumed to be 1 so H0 = sqrt(8*pi/3).

ChaNGa and Gasoline both force G = 1


---------------
Installation
--------------

git clone https://github.com/mtremmel/ICInG.git

cd ICInG

python setup.py install

***or***

python setup.py install --user

***or***

python setup.py install --prefix=<path of your python distribution, e.g. ~/anaconda>


---------------
Requirements
---------------

An updated distribution of python, scipy, numpy (suggestion: anaconda)

pynbody (git@github.com:pynbody/pynbody.git)

tipsy_tools (git@github.com:N-BodyShop/tipsy_tools.git)


---------------
Runnning
---------------

Example (in python environment):

>>import ICInG

>>ICInG.vanilla.genICs()

------------------
Support and contact
------------------
send any concerns/questions/bugs to Michael Tremmel (m.tremmel6@gmail.com, http://www.astro.washington.edu/users/mjt29/)

------------------
Contributing to the Code
-----------------
If you would like to contribute to the code, please create a fork, add your module or edit existing code, and then submit a pull request. (see https://help.github.com/articles/using-pull-requests)

One area of particular interest would be porting the ICs to other formats to be used with other SPH codes.

-----------------
Citation
----------------
I ask that if you use this code to produce data for publications that you acknowledge ICInG and myself (Michael Tremmel) somewhere in your publication.
