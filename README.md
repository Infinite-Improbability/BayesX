# BAYES-X.
*Created by M. Olamaie et al, October 2013*  
*Modified by Ryan Cox, 2022.*

[![Consistency tests](https://github.com/Infinite-Improbability/BayesX/actions/workflows/consistency.yml/badge.svg)](https://github.com/Infinite-Improbability/BayesX/actions/workflows/consistency.yml)
[![Unit tests](https://github.com/Infinite-Improbability/BayesX/actions/workflows/tests.yml/badge.svg)](https://github.com/Infinite-Improbability/BayesX/actions/workflows/tests.yml)

BAYES-X is a Bayesian inference tool for the analysis of X-ray observations of galaxy clusters. For more details please see the BAYES-X paper (arXiv:1310.1885). If you use BAYES-X, please cite the following papers in your publications:

- BAYES-X Paper: arXiv:1310.1885
- MultiNest Papers: arXiv:0704.3704, arXiv:0809.3437 & arXiv:1306.2144


## License

BAYES-X uses MultiNest as the Bayesian inference engine, therefore users are required to agree to the MultiNest licence (given in file called LICENCE in the multinest folder) before using BAYES-X.

## Compiling BAYES-X

Please edit the Makefile in the main directory according to your compilers & library locations. Do 'make all' to compile the code. On successful compilation an executable called 'BayesX' will be created in the bin folder.

### Required Libraries
BAYES-X requires cfitsio, lapack & MultiNest libraries. MultiNest source code is provided in the folder called multinest. Default Makefile should automatically create MultiNest library in the lib folder.

To install dependencies on a Debian based system use
```
sudo apt install liblapack-dev libcfitsio-dev
```

### MPI Support
The code is MPI compatible. In order to disable the MPI parallelization, remove -DMPI compilation flag.


### Compiling with gfortran
You might need to use the `-ffree-line-length-none` flag while compiling the code with gfortran compiler to remove the restriction imposed by gfortran on line length.

## Running BAYES-X
BAYES-X can be run as follows:
```
./bin/BayesX [input-file]
```
in the serial mode, or:

```
mpiexec -np [n] ./bin/BayesX [input-file]
```
in parallel mode. Here n is the number of cores & input-file is a text file with input parameters (read 
the following section). Please see infile64by64by32.inp or infileA262_256by256by433.inp for an example input file.

The datafiles referred to in infile64by64by32.inp and infileA262_256by256by433.inp are provided in the data folder. User can use run BAYES-X with
both input files  to see an example of how it works.

## Input Parameters

In order to run BAYES-X, the user must provide an input file (see e.g. infile64by64by32.inp). This file contains pairs of tags & parameter values where tags identify the inputs & start with the character '#'.
Parameter values corresponding to a given tag are specified in the following line.
Details of these tags are as follows.

### Observation and general configuration
#XrayTelescope:		The first three letters of the X-ray telescope name, e.g. #XrayTelescope = 'CHA' for Chandra, 'XMM' for XMM, etc. Not actually used?  
#nx:			Number of pixels in x direction, e.g. nx = 256  
#ny:			Number of pixels in y direction, e.g. ny = 256  
#xrayNbin:		Number of energy bins, e.g. Nbin = 63  
#xrayNch:		Number of energy channels, e.g. Nch = 63  
#n:			Number of steps for discretising r, e.g. n = 100  
#Aeffave:		Average effective area of the telescope in cm^{2}, e.g. Aeffave = 250 for Chandra. Only used for calculating background.  
#xraycell:		Spatial pixel size in arcsecond, e.g. xraycell = 0.492  
#xrayEmin:		Minimum value of the energy range in keV, e.g. xrayEmin = 0.7  
#xrayEmax:		Maximum value of the energy range in keV, e.g. xrayEmax = 7.0  
#sexpotime:		Source exposure time in second, e.g.sexpotime = 300000  
#bexpotime:		Background exposure time in second, e.g. bexpotime = 300000  
#NHcol:			Hydrogen column density in cm^2, e.g. NHcol = 2.2d20  
#xrayBG_model           Predicted background rate at each pixel in counts cm^-2 arcmin^-2s^-1, e.g. xrayBG_model=8.4d-6

### Radius limits_fraction
#rauto:     Do dynamic limit calculation [T/F]. Overrides following options, default value T  
#rmin:      real, Minimum radius for integration and in sky plane, Mpc  
#rmax:      real, Maximum radius in sky plane, Mpc  
#rlimit:    real, Maximum radius for integrals, Mpc  

### MEKAL
The root for MEKAL data files:  
#filion:		mekal1.dat   
#filrec:		mekal2.dat  
#filkar:		mekal3.dat  
#filcon:		mekal4.dat  
#fillin:		mekal5.dat  
#filcaf:		mekal6.dat  

### Response functions
The root for telescope ARF and RMF files:  
#filARF:		ARF telescope file in .txt format. This is a 1D array of size of number of energy bins (1 times xrayNbin).  
#filRMF:		RMF telescope file in .txt format. This is a 2D array of size of number of energy bins times number of energy channels (xrayNbin times xrayNch).  

### Data file paths
#filBG:			Background file contains the background counts per second per pixel. BAYES-X reads the background counts in a 1-D array of dimension of (nx times ny times xrayNch).
#filevent:		The event file contains the counts per pixel per energy channel. BAYES-X reads the event file in a 1-D array of dimension of (nx times ny times xrayNch) 
#filmask:       The mask file contains 1 where pixels are masked and 0 where they are not. The dimensions match the event file.

### Model selection
#mass_function:	Type of the mass function (with mass_function = 1, 2 & 3 for Evrard, Jenkins & Tinker mass functions respectively). Used with joint M z prior.  
#cluster_model:  Type of the clust model (with cluster_model =1, 2 or 3 for NFW, Einasto or polytropic profiles respectively). The polytropic model is based on Ghirardini2019 [A&A 627, A19 (2019)]

### Priors for the free parameters  
#x_prior:		prior on x coordinate on the sky in arcsec of the cluster center  
#y_prior:		prior on y coordinate on the sky in arcsec of the cluster center  
#m200_prior:		prior on cluster mass within the overdensity radius of r_{200} (M_{200}) in M_{\sun}  
#fgas200_prior:		prior on cluster gas mass fraction within the overdensity radius of r_{200} (f_{g,200})  
#a_GNFW_prior:		prior on slope parameter "a" in GNFW pressure profile  
#b_GNFW_prior:		prior on slope parameter "b" in GNFW pressure profile  
#c_GNFW_prior:		prior on slope parameter "c" in GNFW pressure profile  
#c500_GNFW_prior:	prior on gas concentration parameter, c_{500}  
#alpha_model2_prior:    prior on Einasto shape parameter, alpha  
#gamma0_poly_prior:     prior gamma0, free parameter in polytropic model  
#gammaR_poly_prior:     prior gammaR, free parameter in polytropic model  
#t0_poly_prior:     prior T0, free parameter in polytropic model
#rmin_fraction: prior on minimum integration radius, in units of NFW scale radius r_s. Currently used by model 1 only.
#z_Prior:		prior on cluster redshift  

In the line following the prior tag, users should list 3 values as follows:  
a	b	c  
where a is the prior type, b & c are the parameters specific to the prior type. Allowed prior types & their specific parameters are as follows:

Delta function prior				0   [parameter value]   [parameter value]  
Uniform prior						1   [minimum value]     [maximum value]  
Uniform prior in log				2   [minimum value]     [maximum value]  
Gaussian prior						3   [mean]              [standard deviation]  
Log-normal prior					4   [mean]              [width]  
Joint prior for M_{200} & z using the mass function	8   [minimum value]     [maximum value]  

### MultiNest parameters
#IS:			Do Importance Nested Sampling (INS)? [T/F], default value: F  
#multimodal:		Do mode separation? [T/F], default value: F  
#nlive:			Number of live points, default value: 1000  
#eff:			Target efficiency, default value: 0.8  
#tol:			Tolerance value for convergence, default value: 0.5  
#seed:			Seed for random number generator, negative value would result in seed taken from system clock, default value: -1  
#root:			Root for the output files  
#updint:		No. of iterations after which the output files should be updated, default value: 100  
#maxmodes:		Maximum no. of modes (for memory allocation), default value: 20  
#nCdims:        no. of parameters on which clustering should be performed if mode separation is enabled, default value: 2  

Please refer to the MultiNest README file (provided in multinest folder) for more details on these parameters.

Inputs corresponding to the following tags must be set in the input file otherwise BAYES-X will exit with an error:

#filion  
#filrec  
#filkar  
#filcon  
#fillin  
#filcaf  
#filARF  
#filRMF  
#filBG  
#filevent  
#nx  
#ny  
#xrayNbin  
#xrayNch  
#root  

## Output Files

BAYES-X produces different output files with the information about posterior distributions of cluster parameters, which are in the following order:

### Free parameters
Not all priors are used for every model and may be delta functions so associated posteriors may be absent. This affects the numbering but not the order.
For 

1. cluster x coordinate on the sky in arcsec  
2. cluster y coordinate on the sky in arcsec  
3. cluster redshift, z  
4. cluster mass within the overdensity radius of r_{200} (M_{200}) in M_{\sun}  
5. cluster gas mass fraction within the overdensity radius of r_{200} (f_{g,200})  
6. slope parameter "a" in GNFW pressure profile  
7. slope parameter "b" in GNFW pressure profile  
8. slope parameter "c" in GNFW pressure profile  
9. gas concentration parameter, c_{500}  
10. shape parameter in Einasto profile, alpha  

### Derived parameters
11. angular diameter distance in Mpc
12. scale radius, r_{s} in NFW profile in Mpc or r_2 in Einasto profile in Mpc  
13. the normalisation coefficient in NFW profile, \rho_{s}, in M_{\sun} Mpc-3 or \rho_2 in Einasto profile.  
14. scale radius, rp, in GNFW pressure profile in Mpc  
15. the normalisation coefficient, Pei, in the GNFW pressure profile in keV m^{-3}  
16. cluster overdensity radius of r_{2500} in Mpc  
17. cluster concentration parameter at r_{2500}, c_{2500}  
18. cluster gas mass within within the overdensity radius of r_{2500} (M_{g,2500}) in M_{\sun}  
19. cluster mass within the overdensity radius of r_{2500} (M_{2500}) in M_{\sun}  
20. cluster gas mass fraction within the overdensity radius of r_{2500} (f_{g,2500})  
21. cluster gas temperature at r_{2500} (T_{g,2500}) in keV  
22. cluster electron number density at r_{2500} (n_{e,2500}) in m^{-3}  
23. cluster entropy at r_{2500} (K_{e,2500}) in keVm^2  
24. cluster electron pressure at r_{2500} (P_{e,2500}) in keV m^{-3}  
25. cluster overdensity radius of r_{500} in Mpc  
26. cluster concentration parameter at r_{500}, c_{500}  
27. cluster gas mass within within the overdensity radius of r_{500} (M_{g,500}) in M_{\sun}  
28. cluster mass within the overdensity radius of r_{500} (M_{500}) in M_{\sun}  
29. cluster gas mass fraction within the overdensity radius of r_{500} (f_{g,500})  
30. cluster gas temperature at r_{500} (T_{g,500}) in keV  
31. cluster electron number density at r_{500} (n_{e,500}) in m^{-3}  
32. cluster entropy at r_{500} (K_{e,500}) in keVm^2  
33. cluster electron pressure at r_{500} (P_{e,500}) in keV m^{-3}  
34. cluster overdensity radius of r_{200} in Mpc  
35. cluster concentration parameter at r_{200}, c_{200}  
36. cluster gas mass within within the overdensity radius of r_{200} (M_{g,200}) in M_{\sun}  
37. cluster mass within the overdensity radius of r_{200} (M_{200}) in M_{\sun}  
38. cluster gas mass fraction within the overdensity radius of r_{200} (f_{g,200})  
39. cluster gas temperature at r_{200} (T_{g,200}) in keV  
40. cluster electron number density at r_{200} (n_{e,200}) in m^{-3}  
41. cluster entropy at r_{200} (K_{e,200}) in keVm^2  
42. cluster electron pressure at r_{200} (P_{e,200}) in keV m^{-3}  
43. -146  cluster radius, mass, gas mass, gas mass fraction, gas temperature, electron number density, entropy and pressure at 0.04\timesr_{500} and then at radii of 
(0.05, 0.06,0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)\times r_{500}.

Description of different output files is as follows. Please refer to MultiNest README file (provided 
in the multinest folder) for more details.

[root].txt
Compatible with getdist, this file has 2+nPar columns. Columns have sample probability, -2*loglike,
parameter values. Sample probability is the sample prior mass multiplied by its likelihood & normalized 
by the evidence.

[root]post_separate.dat
This file is only created if multimodal is set to T & contains posterior samples separated by 2 blank lines. 
Format is the same as [root].txt file.

[root]stats.dat
Contains the global log-evidence, its error & local log-evidence with error & parameter means & standard
deviations as well as the  best fit & MAP parameters of each of the mode found with local log-evidence.

[root]post_equal_weights.dat
Contains the equally weighted posterior samples. Columns have parameter values followed by loglike value.

[root]summary.txt
There are nmode+1 (nmode = number of modes) rows in this file. First row has the statistics for the global 
posterior. After the first line there is one row per mode with nPar*4+2 values in each line in this file. 
Each row has the following values in its column mean parameter values, standard deviations of the parameters, 
bestfit (maxlike) parameter values, MAP (maximum-a-posteriori) parameter values, local log-evidence, maximum 
loglike value. If IS = T (i.e. INS being used), firs row has an additional value right at the end, INS 
log-evidence estimate.

### Visualization of output
The [root].txt file created by MultiNest is compatible with the format required by getdist package which is 
part of CosmoMC package. Refer to the following website in order to download or get more information about
getdist:
http://cosmologist.info/cosmomc/readme.html#Analysing

Johannes Buchner's PyMultiNest can also be used on existing MultiNest output to plot & visualize results.
https://github.com/JohannesBuchner/PyMultiNest

## Python Interface

There is an experimental Python wrapper, `pybayesx` to simplify data processing. Documentation can be found [here](https://infinite-improbability.github.io/BayesX/index.html).

## Development Notes

### Profiling
I've used [Score-P](https://www.vi-hps.org/projects/score-p/) and [Cube](https://www.scalasca.org/scalasca/software/cube-4.x/download.html) for profiling. If you have these utilities installed then Bayes-X can be compiled with profiling support using
```
make PREP="scorep"
```
Profiling will then be performed automatically when run and the output can be viewed with Cube.