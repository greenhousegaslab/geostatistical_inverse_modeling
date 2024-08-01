%-----------------------------------------------------------------------------------------------------------------------%
% SCRIPT:  inversion_direct												%
% PURPOSE: Launch a geostatistical inverse model (GIM).									%
%	   This script will estimate the fluxes using the direct equations to solving the GIM.				%
%	   These equqtions can be found in Kitanidis (1996).  								%
%	   I.e., this script does not estimate the fluxes iteratively.							%
% S. Miller, Jan. 8, 2016												%
%	Updated July 24, 2024												%
%															%
%-----------------------------------------------------------------------------------------------------------------------%

%---------------------------%
% Summary of this script    %
%---------------------------%

	% This Matlab script launches the inversion. 
	% Each section of this script should be customized to read in the inputs for your particular inversion.
	% Once you have set up this script, you should be able to launch the inversion without modifying any more scripts.
	% ** I have marked sections with an asterick if that section needs to be customized by the user.
	
%----------------------%
% Set script options   %
%----------------------%

	% Set the following flag to one to estimate the fluxes (i.e., posterior best estiamte)
	estflux   = 1;
	
	% Set the following flag to one to estimate the posterior covariance matrix.
	estuncert = 0;

	
%----------------------------%
% Initial setup parameters   %
%----------------------------%
	
	%---------------------------------------%
	% **Set the name of the output file**   %
	%---------------------------------------%

	% This script will save the estimated fluxes into the following file (csv file):
	outpath = '/home-4/smill191@jhu.edu/work/smiller/GIM_test/';
	outfile = strcat(outpath,'fluxes.csv');
	
	
	%---------------------------------------------%
	% **Set the covariance matrix parameters**    %
	%---------------------------------------------%
	
	% These parameters define the covariance matrices used in the inversion -- R and Q.
	% For more details on R and Q, refer to Anna Michalak's 2004 paper ("A geostatistical approach to surface flux estimation ...")
	
	% We'll refer to the vector of covariance matrix parameters as "theta". Here are the individual elements of theta:
	% theta(1): This value is the square root of the diagonal elements of the R matrix. 
	%	For methane and N2O, theta(1) will typically have units of ppb.
	%	In this setup, I have assumed that R is a diagonal matrix and that the diagonal elements will be constant.
	%	One could certainly re-configure this script to support a more complex formulation of R.
	
	% theta(2): This value is the square root of the diagonal elements of the Q matrix.
	%	This value should have units of micromol m-2 s-1.
	
	% theta(3): This value defines the decorrelation lengthscale in Q. 
	%	This value will typically have units of km.
	
	% theta(4): This value defines the decorrelation timescale in Q.
	%	This value will have units of time. The units will depend on your particular inversion setup.
	%	For example, I typically estimate CH4 fluxes that can vary on a daily time resolution. 
	%	As a result of this setup, I usually define theta(4) to have units of days.
	
	% Set the covariance matrix parameters for the inversion
        theta = [ 2.000 10.046 555.420 9.854 ];
	
	% Display the covariance matrix parameters on screen
	disp('Covariance matrix parameters');
	disp(theta);


%-------------------------------------------------------------%
% Set the grid box areas (for the uncertainty calculations)   %
%-------------------------------------------------------------%

	% Create a vector that contains the area of each grid box in the inverse model.
        % The areas should have units corresponding to the units of the fluxes.
	% The data should be in a column vector with length m. Set a zero in grid boxes that you do not want
        % included in the uncertainty calculations.
	% Put a zero in grid boxes that you do not want included in the uncertainty calculations.       

	 selu = [ FILL IN HERE ];


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

        %! NOTE: See the readme file for more information on different options for handling H.
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

        % Create a vector of ones for the first column of the X matrix
	% [ FILL IN CODE HERE ]

        % APPEND ANY OTHER DESIRED COLUMNS TO X HERE (e.g., an emissions inventory, land surface data, meteorology reanalysis, etc.)


%-------------------------%
% Create the E matrix     %
%-------------------------%

        disp('Create the E and D matrices');

        % The E and D matrices are used to construct the Q covariance matrix.
        % The E matrix describes how covariances in Q decay in space, and the D matrix describes how covariances in Q decay in time.
        % See the following paper by Vineet Yadav for more details on the E and D matrices: http://www.geosci-model-dev.net/6/583/2013/gmd-6-583-2013.html
        % The E matrix has dimensions (r x r)
        % The D matrix has dimensions (q x q)


        %----------------------------------------------%
        % **Read in the geographic distance matrix**   %
        %----------------------------------------------%

        % This matrix (dimensions r x r) should define the distance between each emissions grid box
        % This matrix is symmetric with zeros on the diagonals.
        % This matrix should typically has units of km. Whatever units you choose should match the units of theta(3)

        deltamat = [ READ IN THE MATRIX HERE ];


        %------------%
        % Create E   %
        %------------%

        % No need to edit this section
        % Note: I use a spherical covariance model because of its sparse properties.

        % Multiply the distance matrix by the covariance model
        % Exponential covariance model
        % decaymat = exp(-1 .* deltamat ./ theta(3));

        % Spherical covariance model
        E                      = 1 - 1.5 .* (deltamat ./theta(3))  + 0.5 .* (deltamat.^3 ./ theta(3).^3);
        E(deltamat > theta(3)) = 0;
	

%----------------------------%
% Create the D matrix        %
%----------------------------%

        % The D matrix describes the time between each time period in the inversion
        % It will have dimensions (q x q)
        % The D matrix is usually easier to set up than the E matrix, and I have included sample code below

        %-----------------------------------%
        % **Create time distance matrix**   %
        %-----------------------------------%

        ntimes = [ FILL IN HERE WITH THE NUMBER OF TIME PERIODS IN YOUR INVERSION ];
        days = 1:ntimes;
        days = days';
        days = days * ones(1,length(days));
        days = abs(days - days');


        %------------%
        % Create D   %
        %------------%

        % No need to edit this section
        % Exponential model
        % D      = exp(-1 .* days ./ theta(4));

        % Spherical covariance model
        D = 1 - 1.5 .* (days   ./theta(4))  + 0.5 .* (days.^3   ./ theta(4).^3);
        D(days   > theta(4)) = 0;


%----------------------------%
% Create the sigmaQ vector   %
%----------------------------%

        % This section defines the diagonal elements of the Q matrix.
        % There is no need to edit this section if you want Q to have constant diagonal elements.
        % If you want Q to have varying diagonal elements, then you will need to edit this section.

        ntimes = size(D,1);
        sigmaQ = theta(2) .* ones(ntimes,1);

        % Multiply the sigmaQ vector by D
        % Note: in the line of code below, I assume that sigmaQ may be different at different times.
        % However, I asssume that sigmaQ does not change spatially.
        % One could construct an inverse model where sigmaQ changes spatially, but you would need to edit
        % the code below accoridngly.
        D = (sigmaQ*sigmaQ') .* D;


%--------------------------------%
% Create Q using D and E         %
%--------------------------------%

	Q = kronMat(D,E);


%------------------------%
% Create the R matrix    %
%------------------------%
		
	disp('Create R');
	tic; 
	
	% No need to edit this section, unless you want a more complicated setup for R.
	R = (theta(1).^2) .* ones(n,1);
	R = spdiags(R,0,n,n);

	disp(toc);

        % Define the sizes of the different matrices
        p   = size(X,2);
        n   = length(Z);
        m1  = size(E,1);
        ntimes = size(D,1);
	m   = ntimes .* m1;


%---------------------------------------------------%
% Pass the columns of X through the forward model   %
%---------------------------------------------------%

        disp('Calculate HX');
	H = matvecH(ntimes,Hpath);
	HX = H*X;


%---------------------------%
% Compute the psi matrix    %
%---------------------------%

        disp('Launch the script to calculate HQH + R');
        psi = HQHR(R, D, E, Hpath);


%-------------------------------%
% Estimate the fluxes		%
%-------------------------------%

	if estflux == 1;

	disp('Compute the left hand side of the inverse modeling equations');
	LHS = [[psi HX] ; [HX' zeros(p,p)]];
	if estuncert==0; clear psi HX; end;

	% Compute the inverse modeling weights
	disp('Compute the weights');
	RHS = [ Z ; zeros(p,1)];
	weights = LHS \ RHS;

	% Separate out the individual components of the 'weights' vector 
	eta   = weights(1:(size(weights,1)-p));
	betas = weights((size(weights,1)-p+1):size(weights,1));;	

	% Calculate H * eta
	Heta = H' * eta;

	% Calculate Q * H^T * eta
	disp('Calculate Q * H^T * eta');
	QHeta = Q * Heta; 

	% Create the estimate of the fluxes
	disp('Estimate the fluxes');
	shat = X*betas + QHeta;


%------------------------------%
% Write the outputs to file    %
%------------------------------%

        disp('Writing outputs to file');
        dlmwrite(outfile,full(shat),',');

        disp('Outputs written to file');
        disp(outfile)

        end; % End of estflux if statement


%-------------------------------%
% Calculate the uncertainties   %
%-------------------------------%

	if estuncert == 1;

        uncert  = uncert_direct(psi,HX,X,Q,H,selu);

	end; % End of estuncert if statement


%-----------------------------------------------------------------------------------------------------------------------%
% END OF SCRIPT														%
%-----------------------------------------------------------------------------------------------------------------------%
