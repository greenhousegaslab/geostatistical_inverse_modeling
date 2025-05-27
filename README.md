# geostatistical_inverse_modeling
This repository contains Matlab code for conducting geostatistical inverse modeling using large atmospheric datasets.

Authors: Scot Miller and Arvind Saibaba. (Note: we want this code repository to be a collaborative effort, and the more authors/contributors, the better!)

* This code has been updated to include new class definitions for Q and H. These new class definitions substantially simplify the code.
(Class definitions are adopted from Julianne Chung and Taewon Cho.)

----------------------------------------------------------
Introduction and overview
----------------------------------------------------------

This code repository contains several Matlab scripts for running a geostatistical inverse model (GIM). We chose Matlab because it is a scripting language (versus a compiled language) and because matrix algrebra operations are generally faster in Matlab than other scripting languages like R or Python. This repository includes several different algorithms for solving a GIM, and these algorithms should all converge on the same answer. The choice of algorithm will likely depend upon the size and characteristics of your inverse problem (see below for more detail).

Most aspects of these scripts do not need to be customized or edited by the user. However, a few aspects of the scripts need to be customized for your particular inverse modeling setup, and those sections of the scripts are clearly marked. Most importantly, you will likely need to customize these scripts to match your specific atmospheric transport model. Most of the scripts here can be run with both particle-following or Lagrangian models (e.g., STILT for FLEXPART) and gridded or Eulerian models (e.g., GEOS-Chem or TM5; provided that a model adjoint is available). However, you will likely need to either modify the GIM code or modify your atmospheric model output to make the two compatible.

We have included a case study that can be run out-of-the-box (described in more detail below). The purpose of this case study is to show how the code works and how it can be paired with different inverse modeling inputs.

----------------------------------------------------------
Fair use
----------------------------------------------------------

Please send us an email (scot.m.miller [at] gmail.com and asaibab [at] ncsu.edu) if you are thinking of using these scripts in your research. We would love to hear more about what you are working on and if we can help in any way. 

Also, please make sure to cite the companion article in Geoscientific Model Development if you do use these scripts in your work.


----------------------------------------------------------
A note on the structure of the scripts
----------------------------------------------------------

These scripts provide four different ways to run a geostatistical inverse model (GIM).
(1) Solve the GIM equations directly: The direct or analytical solution involves solving a system of linear equations.
This approach works (a) if the number of observations is less than ~50,000 and (b) if an explicit H matrix is readily available. To solve the GIM in this way, modify and run the script inversion_direct.m. This script also provides an option for computing the posterior uncertainties directly.

(2) Estimate the fluxes using the minimum residual approach: In this approach, the fluxes are estimated by iterating toward the solution using a minimum residual method. This approach works well if (a) the number of observations is large (e.g., > 50,000) and/or (b) if an explicit H matrix is not readily available (i.e., if one is using an adjoint model that does not produce an explicit H matrix). Note that this approach is not recommended if the fluxes have known bounds (e.g., if the fluxes cannot be negative). To solve the GIM in this way, modify and run the script inversion_dual_eq.m. This script also provides an option for approximating the uncertainties using a reduced rank approach.

(3) Estimate the fluxes using a limited memory Broyden–Fletcher–Goldfarb–Shanno (L-BFGS) algorithm with a variable transformation: In this approach, the fluxes are estimated by iterating toward the solution using a quasi-Newton algorithm known as L-BFGS.  This approach works well if (a) the number of observations is large (e.g., > 50,000) and/or (b) if an explicit H matrix is not readily available (i.e., if one is using an adjoint model that does not produce an explicit H matrix). This approach should exhibit similar computational performance to the minimum residual approach (2). However, this script does require inverting the covariance matrix R and the components of the Q covariance matrix (D and E). If these matrix inverses are burdensome to compute, the minimum residual approach is likely the better option. The L-BFGS is used in many 4-D variational (4D-Var) inverse models; these scripts may therefore be easier to integrate into existing 4D-Var code. Note that this approach is not recommended if the fluxes have known bounds (e.g., if the fluxes cannot be negative). To solve the GIM in this way, modify and run the script inversion_LBFGS.m.

(4) Estimate the fluxes using the L-BFGS algorithm with bounds: This approach is similar to L-BFGS but will additionally enforce bounds on the inverse problem. For example, this option may be a good choice for methane or N2O, gases that do not have large negative fluxes. This script also iterates toward the solution. However, it converges much more slowly than either options (2) or (3). For example, this script may require ten times as many iterations to converge on a solution, and it is therefore not recommended for very large inverse problems. To solve the GIM in this way, modify and run the script inversion_LBFGS.m. One would need to remove the variable transformation from this script, swap out calls to the functioncost_gradient_fun_transform.m with calls to the function cost_gradient_fun.m, and swap out calls to the L-BFGS algorithm with calls to an L-BFGS-B algorithm (the bounded version). 


----------------------------------------------------------
Which set of scripts should I use?
----------------------------------------------------------

Q: Are there bounds on the fluxes? For example, do the fluxes need to be non-negative? <br>
Use the L-BFGS-B (option 4). One could also use option (1) with a variable transformation to enforce non-negativity (see Snodgrass and Kitanidis 1997 [https://doi.org/10.1029/96WR03753] for one example).

Q: Is the number of observations less than ~50,000, and is an explicit H matrix available? <br>
Use the direct solution (option 1). Sometimes, inverse problems with between 25,000 and 50,000 observations are too large for the direct solution (depending upon the problem and the computer system involved), and one may need to use option 2 or 3.

Q: What if my atmospheric model does not produce an explicit H matrix? <br>
Use option (2) or (3) (or option 4 for a problem with bounds). Note that you will need both a forward model and an adjoint model to use any of the inverse modeling approaches outlined above.


----------------------------------------------------------
A note on the H matrix
----------------------------------------------------------

Some available atmospheric models will produce an explicit H matrix. These models include STILT (the Stochastic, Time-Inverted, Lagrangian Transport Model) and FLEXPART (the FLEXible PARTicle dispersion model). Other atmospheric models do not produce an explicit H matrix (e.g., TM5 or GEOS-Chem). Instead, these models can pass a vector through the forward atmospheric model, and in some cases, pass a vector through the adjoint model. Specifically, these models will calculate the product H or H^T and a vector but will not formulate H or H^T explicitly.

You will need to modify the scripts here depending upon which type of model that you use. Most of these scripts can be used with either type of model.

The scripts here have been written for a model that produces an explicit H matrix, but these scripts can be modified to account for other types of atmospheric models. 

The scripts included here are written assuming the following format for the H matrix:
- An explicit H matrix is available.
- The H matrix has been broken into smaller pieces or strips, and each strip has been written to a different file. Each strip should have n rows (corresponding to the total number of observations in the inverse model), and each strip should have columns corresponding to a single time period of fluxes to be estimated. 
- These H "strips" are saved in Matlab files (with the ".mat" extension). The files should have the following naming convention: "H_1.mat" where 1 indicates the time period that this particular file corresponds to. In this case, "H_1.mat" would refer to the first time period and "H_2.mat" would correspond to the second time period of the inverse model, etc. For example, if one were estimating 3-hourly CO<sub>2</sub> fluxes, H_1.mat would correspond to the first three-hourly time period, and H_2.mat to the second three-hourly time period.
- Again, this setup can be modified for atmospheric models that do not produce an explicit H matrix or for inverse problems where the H matrix is stored in a different format. I have set up the scripts this way solely out of convenience; the setup described above lends itself to the case study included here using the STILT model. I have flagged all of the sections of code that would need to be modified to run the GIM with different types of atmospheric models that do not fit the H matrix format described above. Look for sections of code denoted with the comment "%! Note: Edit this section to match the actual format of H in your problem."

----------------------------------------------------------
Computational considerations
----------------------------------------------------------

The computational approaches implemented in this code are designed for large inverse problems, but it is nevertheless important to keep computational considerations in mind when adapting the code for a specific inverse problem. We discuss several of these considerations below:

(1) The number of iterations required by the iterative solver to estimate the fluxes can be an important limiting factor when using certain types of adjoint atmospheric models but may not be a limiting factor when using other types of atmospheric models. For trajectory models like the Stochastic Time-Inverted Lagrangian Transport (STILT) model, <b>H</b> is formulated explicitly and can be read in directly. In that case, the computing resources required to run numerous STILT trajectories, not the number of iterations required by the solver, is likely to be the computational bottleneck. By contrast, the number of iterations required to converge on a solution is likely to be the bottleneck for gridded chemical transport models like GEOS-Chem or TM5. These models do not produce an explicit <b>H</b> and <b>H</b><sup>T</sup> matrices, and one must instead run the forward and adjoint models once per iteration of the solver. These calculations often become time-intensive when numerous iterations are required to converge on a flux estimate. Furthermore, some adjoint models (i.e., GEOS-Chem) cannot be run in parallel for greenhouse gas applications, though we expect that these capabilities will change in the future with the development of an adjoint for models like GEOS-Chem-High Performance (GC-HP).

(2) The matrices <b>D</b> and <b>E</b> (Eqs. 9-10 in the GMD manuscript) are usually straightforward to store in memory and/or invert given the dimensions of most atmospheric inverse models to date. However we anticipate that this will change in the future as atmospheric models like GEOS-Chem have better parallel computing capabilities and can be run at higher spatial resolution. In those cases, it may be important to structure <b>D</b> and <b>E</b> as hierarchical matrices or structured matrices (e.g., Toeplitz matrices) to avoid problems with storing these matrices in memory or inverting these matrices. Hierarchical and structured matrices are not implemented in this code repository, but Sect. 3 of the GMD manuscript provides several references on this topic.

(3) The choice of covariance function can have a large impact on the wall clock time and memory required for matrix calculations using <b>D</b> and <b>E</b>. An exponential covariance model is very common in existing GIM studies in hydrology and atmospheric science. For large inverse problems, this choice may not be practical; an exponential model will never decay to zero. As a result, <b>D</b> and <b>E</b> will never be sparse matrices. By contrast, other covariance models, like a spherical model, do decay to zero, and <b>D</b> and <b>E</b> can be formulated as memory-saving sparse matrices.

(4) The code here can be re-written for other languages if a different language is more convenient than Matlab. We recommend that users exercise caution if doing so because different commonly-used languages can exhibit very different performance. For example, we found that R is far slower than Matlab at linear algebra and often requires more memory than Matlab for the same matrix inversion.


----------------------------------------------------------
List of scripts in the repository
----------------------------------------------------------

* OPTION 1: Solve the GIM equations directly

SCRIPT: inversion_direct.m <br>
PURPOSE: Estimate the fluxes using the direct solution (i.e., by solving a system of linear equations). <br>
CALLED BY: None. <br>
CALLS: HQHR.m, uncert_direct.m

SCRIPT: HQHR.m <br>
PURPOSE: Calculate (H*Q*H' + R). This matrix is often referred to with the Greek letter Psi. <br>
CALLED BY: inversion_direct.m <br>
CALLS: None.

SCRIPT: uncert_direct.m <br>
PURPOSE: Calculate the posterior uncertainties using direct calculations. This script will provide an exact answer
	 and will not use any approximations to improve computational tractability. <br>
CALLED BY: inversion_direct.m <br>
CALLS: None.


* OPTION 2: Minimum residual approach

SCRIPT: inversion_dual_eq.m <br>
PURPOSE: Launch the geostatistical inverse model using the minimum residual approach to estimating the fluxes.
	 This script also provides an option for estimating the uncertainties using a reduced rank approach. <br>
CALLED BY: None. <br>
CALLS: Ax.m, weights_to_fluxes.m, uncert_est.m

SCRIPT: Ax.m <br>
PURPOSE: Calculate the left-hand side of the GIM equations used in the minimum residual approach. <br>
CALLED BY: inversion_dual_eq.m <br>
CALLS: None.

SCRIPT: weights_to_fluxes.m <br>
PURPOSE: Convert the weights estimated using the minimum residual approach into estimated fluxes. <br>
CALLED BY: inversion_dual_eq.m <br>
CALLS: None.

SCRIPT: uncert_est.m <br>
PURPOSE: Estimate uncertainties using the reduced rank approach. This approach estimates approximate
	 uncertainties. It is computationally tractable even for very large inverse problems. <br>
CALLED BY: inversion_dual_eq.m <br>
CALLS: QBABQ.m, randeigdecomp.m, eigfun.m

SCRIPT: QBABQ.m <br>
PURPOSE: Calculate Q^(1/2)*B*A*B'*Q^(1/2)*B where B is some vector or matrix. <br>
CALLED BY: uncert_est.m <br>
CALLS: None.

SCRIPT: eigfun.m <br>
PURPOSE:  Function for computing Ax where A is the matrix we want eigenvectors/values for.  <br>
CALLED BY: uncert_est.m, randeigdecomp.m <br>
CALLS: None.

SCRIPT: randeigdecomp.m <br>
PURPOSE: Compute the eigenvectors and eigenvalues of matrix A using a randomized approach. <br>
CALLED BY: uncert_est.m <br>
CALLS: eigfun.m


* OPTION 3-4: L-BFGS and L-BFGS-B

SCRIPT: inversion_LBFGS.m <br>
PURPOSE: Estimate fluxes using a geostatistical inverse model, implemented using 
	the L-BFGS algorithm and a variable transformation. <br>
CALLED BY: None. <br>
CALLS: cost_gradient_fun_transform.m

SCRIPT: inversion_LBFGS_bounded.m <br>
PURPOSE: Estimate fluxes using a geostatistical inverse model, implemented using 
	the L-BFGS algorithm. This script will also enforce bounds on the solution
 	I.e., it will ensure that the estimated emissions are not negative. <br>
CALLED BY: None. <br>
CALLS: cost_gradient_fun.m

SCRIPT: cost_gradient_fun_transform.m <br>
PURPOSE: Calculate the cost function and gradient for the geostatistical inverse model. <br>
CALLED BY: inversion_LBFGS.m <br>
CALLS: None.

SCRIPT: cost_gradient_fun.m <br>
PURPOSE: Calculate the cost function and gradient function without a variable transformation. This function is an alternative to the script cost_gradient_fun_transform.m. Specifically, this script is useful for inverse problems with bounds (e.g., non-negativity). <br>
CALLED BY: inversion_LBFGS_bounded.m <br>
CALLS: None.


----------------------------------------------------------
Scripts customized to run the included case study
----------------------------------------------------------

The case study here is identical to the first case study in Miller et al. (2019). As a result, the fluxes estimated by these scripts can be compared directly with the first case study in Miller et al. (2019). All of the input files for the case study are available on Zenodo at www.dx.doi.org/10.5281/zenodo.3241467.

In this case study, we estimate 6 weeks of CO<sub>2</sub> fluxes (July through mid-August 2015) across North America using observations from the OCO-2 satellite. Specifically, these scripts will estimate CO<sub>2</sub> fluxes at a 3-hourly time resolution and 1 degree by 1 degree spatial resolution across the US. The observations used in this case study are synthetic; they were generated using fluxes from NOAA's CarbonTracker product and using the Weather Research and Forecasting (WRF) model paired with the Stochastic Time-Inverted Lagrangian Transport (STILT) model. We further added noise to the synthetic observations to give them qualities more similar to real-world observations.

We use an explicit H matrix for this case study; these inputs were generated using the WRF-STILT model as part of NOAA's CarbonTracker-Lagrange project (https://www.esrl.noaa.gov/gmd/ccgg/carbontracker-lagrange/). 

Here is a summary of the different scripts that can be used to run this case study. To run each of these scripts, one only needs to provide a path to the model and data input files; all of the other inputs should be customized to the case study. These scripts will write the estimated fluxes (averaged over the entire 6-week time period) as a netcdf file. This file should contain all of the information needed to plot the estimated fluxes. 
 
- casestudy_inversion_direct.m: Estimate fluxes and uncertainties using the direct approach.

- casestudy_inversion_dual_eq.m: Estimate fluxes and uncertainties using the minimum residual method. Also estimate uncertainties using a reduced rank approach.

- casestudy_inversion_LBFGS.m: Estimate fluxes and uncertainties using the L-BFGS algorithm.

----------------------------------------------------------
References
----------------------------------------------------------

Miller, S. M., Saibaba, A. K., Trudeau, M. E., Mountain, M. E., and Andrews, A. E.: Geostatistical inverse modeling with very large datasets: an example from the Orbiting Carbon Observatory 2 (OCO-2) satellite, Geosci. Model Dev., 13, 1771–1785, https://doi.org/10.5194/gmd-13-1771-2020, 2020.

----------------------------------------------------------
Contact information
----------------------------------------------------------

Scot Miller <br>
Johns Hopkins University <br>
Department of Environmental Health and Engineering <br>
Email: smill191 [at] jhu.edu OR scot.m.miller [at] gmail.com <br>
https://scotmmiller.github.io/ 

Arvind Saibaba <br>
North Carolina State University <br>
Depatment of Mathematics <br>
Email: asaibab [at] ncsu.edu <br>
https://asaibab.math.ncsu.edu/



