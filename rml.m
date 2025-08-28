%---------------------------------------------------------------------------------------------------------------------------------------%
% FUNCTION: rml.m															%
% PURPOSE: Estimate the covariance matrix parameters for a GIM using the restricted maximum likelihood.					%
% S. Miller, July 27, 2015														%
% 	Updated Aug. 28, 2025														%
%																	%
%---------------------------------------------------------------------------------------------------------------------------------------%

%-----------%
% Notes:    %
%-----------%

        % We'll refer to the vector of covariance matrix parameters as "theta". Here are the individual elements of theta:
        % theta(1): This value is the square root of the diagonal elements of the R matrix. 
        %       For methane and N2O, theta(1) will typically have units of ppb.
        %       In this setup, I have assumed that R is a diagonal matrix and that the diagonal elements will be constant.
        %       One could certainly re-configure this script to support a more complex formulation of R.

        % theta(2): This value is the square root of the diagonal elements of the Q matrix.
        %       This value should have units of micromol m-2 s-1.

        % theta(3): This value defines the decorrelation lengthscale in Q. 
        %       This value will typically have units of km.

        % theta(4): This value defines the decorrelation timescale in Q.
        %       This value will have units of time. The units will depend on your particular inversion setup.
        %       For example, I typically estimate CH4 fluxes that can vary on a daily time resolution. 
        %       As a result of this setup, I usually define theta(4) to have units of days.

	% This script is written to estimate a single diagonal value for R, a single diagonal element for Q, a single
	% 	decorrelation length, and a single decorrelation time (as outlined above). R and Q can have any structure. 
	%	If you want to structure R and Q differently from the simple structure above, you would need to modify this script
	% 	accordingly.

	% This script uses a spherical covariance model for the spatial and temporal covariances in Q. 
	%	This model is convenient because of its sparse properties.
	%	If you decide to use a difference covariance model, then modify this script accordingly.
	%	An exponential model is also a common choice for inverse modeling.

	% In some cases, RML isn't sensitive to certain covariance parameters. I.e., in some cases, the cost function will
	%	change very little in response to large changes in one of the parameters, and the procedure here will not converge
	%	on a stable value for that parameter.  If this happens, it usually indicates
	%	that the observations (Z) are not very sensitive to the choice of covariance parameter. In this case, it often
	%	makes more sense to estimate that value based upon expert judgement and optimize the remaining parameters using RML. 

	% The function below will compute the RML cost function. The best covariance matrix parameters are those that minimize
	%	the cost function below. To minimize this function, I use the following line of code in Matlab:
	%	theta = fminsearch( @(theta) reml(theta, Z, X, H, obsites, deltamat, months), theta0);


%----------------------------------------------------------------%
% ** Set an initial guess for the covariance matrix parameters** %
%----------------------------------------------------------------%

	% Note: I would modify this initial guess for your particular setup. 
	% (The initial guess here may not be a good guess for your particular problem.)
	theta0  = [ 16 0.128 189.7 105 ];


%---------------------------------%
% **Read in the observations**    %
%---------------------------------%

        disp('Read in the observation vector');

        % Read in the vector of observations here.
        % Subtract off the boundary condition before reading in the observations here.  
        % We'll refer to this vector as "Z".
        % This vector has length (n x 1) where n are the number of observations

        Z = [ READ IN OBSERVATIONS HERE ];
	n = length(Z);



%-------------------------------%
% **Read in the footprints**    %
%-------------------------------%

        disp('Read in the footprints');

        % Read in the matrix of footprints here (H).
        % The footprint matrix should have n rows, and the rows should correspond to the observations in Z.
        % The footprint matrix should have (q x r) columns.
        % q = the number of time periods in the inversion
        % r = the number of flux grid boxes in the geographic domain
        % The columns of the H matrix should be ordered as follows:
        % The first r columns should be from time period q(1). The next r columns should be from time period q(2), etc.

        % If the H matrix is large, you may want to use the function sparse() to save memory (H is almost always a sparse matrix.).

        Hpath = [ SET THE PATH TO THE FOOTPRINT MATRIX HERE ];


%----------------------------%
% **Create the X matrix**    %
%----------------------------%

        disp('Create the X matrix');

        % The X matrix contains spatial and temporal patterns that will be used in the prior model of the inversion.
        % In a GIM, we usually refer to this prior model as the "deterministic model"
        % The X matrix should have (q x r) rows. These rows should correspond to the columns of H.
        % The X matrix will have p columns. Each column contains land surface, environmental, or emissions information that we want to use in the deterministic model.

        % A Bayesian inversion or GIM always assumes that the prior or deterministic model is unbiased relative to the observations. If the prior model is biased relative
        % to the observations, then we have violated the statistical assumptions of the inversion.
        % Thankfully, the GIM has a built-in mechanism to ensure that the deterministic model is not biased relative to the observations. 
        % The GIM will automatically scale all of the components in X to match the actual observations. As a result, the actual magnitude of X is not important in the inversion.

        % Options for X:
        % OPTION 1: In the simplest case, X could be a vector of ones. In this case, X will have p=1 columns. This setup is similar to what Kim Mueller had done in her 2008 CO2 inversion paper. 
        %       If you set up the X matrix in this way, you are making as few assumptions about the fluxes as possible. Any spatial/temporal patterns in the fluxes will come directly
        %       from the observations and not from an emissions inventory or bottom-up emissions model.
        % OPTION 2: One could include an existing emissions inventory in the X matrix (e.g., the EDGAR inventory). In this setup, X would have p=2 columns. The first column would be a 
        %       vector of ones, and the 2nd column would be the chosen emissions inventory. 
        %       In a GIM, it is customary to include a vector of ones in X, irrespective of 
        %       what else gets included in X. This vector of ones is analagous to the intercept term in a regression. In a GIM, this vector of ones will account for any patterns in the fluxes
        %       that aren't captured by other columns of the X matrix. 
        % OPTION 3: X could have a more complex format with land surface data, meteorological reanalysis, etc. Sharon Gourdji's 2008 and 2012 papers provide good examples. 

        % ADD CODE HERE TO CREATE X

%-------------------------------------------------------%
% Create the spatial and temporal distance matrices     %
%-------------------------------------------------------%

        disp('Create the spatial and temporal distance matrices');

        % These matrices are used to construct the Q covariance matrix. 


        %----------------------------------------------%
        % **Read in the geographic distance matrix**   %
        %----------------------------------------------%

        % This matrix (dimensions r x r) should define the distance between each emissions grid box
        % This matrix is symmetric with zeros on the diagonals. 
        % This matrix should typically has units of km. Whatever units you choose should match the units of theta(3)

        deltamat = [ READ IN THE MATRIX HERE ];


        %-----------------------------------%
        % **Create time distance matrix**   %
        %-----------------------------------%

	% Note that I've named the variable here "days", but it's ok to use a different time unit -- so long as "times" and 
	% the theta values have the same units.
        ntimes = [ FILL IN HERE WITH THE NUMBER OF TIME PERIODS IN YOUR INVERSION ];
        days = 1:ntimes;
        days = days';
        days = days * ones(1,length(days));
        days = abs(days - days');


%-----------------------------------%
% Minimize the RML cost function    %
%-----------------------------------%

	% Call the "metvecH" function to create an H matrix object.
        ntimes = size(D,1);
	H = matvecH(ntimes,Hpath);

	% Call "fminsearch" to minimize the RML cost function.
	theta = fminsearch( @(theta) rml_costfun(theta, Z, X, H, deltamat, days), theta0);

	disp('Final estimated covariance parameters');
	disp(num2str(theta));


%---------------------------------------------------------------------------------------------------------------------------------------%
% END OF SCRIPT																%
%---------------------------------------------------------------------------------------------------------------------------------------%
