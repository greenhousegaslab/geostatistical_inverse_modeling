%-----------------------------------------------------------------------------------------------------------------------%
% SCRIPT:  inversion_dual_eq                                                                                            %
% PURPOSE: Launch a geostatistical inverse model (GIM).                                                                 %
%          Estimate fluxes using a minimum residual approach to iterate toward the solution.                            %
% S. Miller, Jan. 8, 2016                                                                                               %
%                                                                                                                       %
%-----------------------------------------------------------------------------------------------------------------------%


%---------------------------------%
% REQUIRED USER-DEFINED OPTIONS   %
%---------------------------------%

	% Set what you want to calculate (e.g., best estimate, uncertainties). See description of options below.
	scriptoption = 2;

        % scriptoption = 1: estimate fluxes with min residual method (no preconditioner)
        % scriptoption = 2: estimate uncertainties with reduced rank approach

        % Set the path to all of the files for the case study
        % The estimated fluxes will also be saved to this folder.
        inpath = '/home-4/smill191@jhu.edu/work/smiller/GIM_test/case_study/';

        % Set the name of the output file
        outfile = strcat(inpath,'fluxes_minres.nc');

        % If scriptoption==1: Set the number of iterations for the minimum residual algorithm
	maxit = 50;

	% If scriptoption==2: Set the maximum number of eigenvectors and values to use
	eigmax = 100;

	% If scriptoption==2: Create a file that contains the area of each grid box in the inverse model.
	% The data should be in a column vector with length m. Set a zero in grid boxes that you do not want
	% included in the uncertainty calculations.
        % Read in a vector that contains grid box areas for the continental US (all other boxes are set to zero in the vector).
        load(strcat(inpath,'areas_us.mat'));

        % The estimated fluxes have units of micromol m-2 s-1.
        % We'll want to convert the units from micromol m-2 s-1 to micromol m-2.
        % In this case study, the inverse model will estimate fluxes for 3-hourly time periods,
        % so we'll convert from (1/s) to (1/3 hours)
        selu = sel.*60.*60.*3;


%-----------------------------------%
% READ IN INVERSE MODELING INPUTS   %
%-----------------------------------%

%---------------------------------------------%
% **Set the covariance matrix parameters**    %
%---------------------------------------------%
	
        % These parameters define the covariance matrices used in the inversion -- R and Q.
        % For more details on R and Q, refer to Anna Michalak's 2004 paper ("A geostatistical approach to surface flux estimation ...")

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

        % Set the covariance matrix parameters for the inversion
        theta = [ 2.000 10.046 555.420 9.854 ];

        % Display the covariance matrix parameters on screen
        disp('Covariance matrix parameters');
        disp(theta);


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

        Hpath = inpath;


%---------------------------------%
% **Read in the observations**    %
%---------------------------------%

        disp('Read in the observation vector');

        % Read in the vector of observations here.
        % Subtract off the boundary condition before reading in the observations here.
        % We'll refer to this vector as "Z".
        % This vector has length (n x 1) where n are the number of observations

        load(strcat(inpath,'Z.mat'));
        n = length(Z);


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

        % In the case study here, X will have 8 columns. Each column corresponds to a different 3-hourly time period of the day (e.g., column 1 corresponds
        % to fluxes from 0 to 3 AM UTC, and the second column corresponds to fluxes from 3 AM to 6 AM UTC).
        X = [] ;
        for i = 1:8
                for j = 1:3222*41*8
                        if rem((fix((j-1)/3222) + 1 - i), 8) == 0
                        X(j, i) = 1 ;
                        else
                        X(j, i) = 0 ;
                        end
                end
        end

        X = sparse(X);


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

        load(strcat(inpath,'distmat.mat'));


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

        % Take the inverse of the E matrix
        Einv  = inv(E);


%----------------------------%
% Create the D matrix        %
%----------------------------%

        % The D matrix describes the time between each time period in the inversion
        % It will have dimensions (q x q)
        % The D matrix is usually easier to set up than the E matrix, and I have included sample code below

        %-----------------------------------%
        % **Create time distance matrix**   %
        %-----------------------------------%

        ntimes = 328;
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
        Dinv = inv(D);


%----------------------------------------------%
% Create Q and Q^(-1) using the kronMat class  %
%----------------------------------------------%

	Q    = kronMat(D,E);
	Qinv = kronMat(Dinv,Einv);


%------------------------%
% Create the R matrix    %
%------------------------%
		
        disp('Create R');

        % No need to edit this section, unless you want a more complicated setup for R.
        R = (theta(1).^2) .* ones(n,1);
        R = spdiags(R,0,n,n);


%-------------------------------%
% Precompute several objects    %
%-------------------------------%

	% Details on the Matlab minimum residual function:
	% x = minres(A,b) attempts to find a minimum norm residual solution x to the system of linear equations A*x=b
	% A can also be a function handle that returns A*x

        % Define the sizes of the different matrices
        p      = size(X,2);
        n      = length(Z);
        m1     = size(E,1);
	ntimes = size(D,1);
	m      = m1 .* ntimes;

	disp('Problem size');
	disp('No of unknown fluxes');
	disp(num2str(m));
	disp('No of fluxes per time increment');
	disp(num2str(m1));
	disp('Number of obs');
	disp(num2str(n));

	%-------------------%
	% Calculate b       %
	%-------------------%

        b = [ Z ; zeros(p,1)];


	%----------------%
	% Calculate HX   %
	%----------------%

        disp('Calculate HX');
	H = matvecH(ntimes,Hpath);
	HX = H*X;


%-----------------------------------%
% OPTION 1: Estimate the fluxes     %
%-----------------------------------%

if scriptoption == 1;

        disp('Estimate weights using minimum residual method');
        tol = 1e-30; % Tolerance
        f1 = @(weights) Ax(R, X, HX, Q, H, weights);
        [weights1, flag,relres,iter,resvec] = minres(f1,b,tol,maxit);


        % Estimate the fluxes using the kriging weights
        shat = weights_to_fluxes(weights1,H,Q,X);


%------------------------------%
% Write the outputs to file    %
%------------------------------%

        disp('Reformat outputs and write to file');

        % Convert the estimated fluxes from a vector to a grid
        % Read in the land mask
        load(strcat(inpath,'land_mask.mat'));

        % Only keep fluxes from the month of July and remove fluxes from June
        % The flux estimate begins on 2015-6-21, and we only want to keep fluxes beginning on July 1, 2015
        % (We included the last week of June in the inverse model primarily to ensure that there were not any edge effects
        % at the beginning of the inverse modeling time preiod).
        % Select out the fluxes from July 1 onward
        sel  = (length(land_mask).*8.*10 + 1):length(shat);
        shat = shat(sel);

        % Average the estimated fluxes across all time periods of the inverse model
        shat = reshape(shat,m1,length(shat)./m1);
        shat = mean(shat,2);

        % Use the land mask to put the flux estimate on a latitude-longitude grid
        lat               = 10.5:79.5;
        lon               = -179.5:-10.5;
        fluxes            = zeros(length(lon),length(lat));
        fluxes(land_mask) = shat;

        % Write the outputs to netcdf file
        nccreate(outfile,'fluxes','Dimensions', {'nlon',size(fluxes,1),'nlat',size(fluxes,2)},'FillValue','disable','Datatype','double');
        nccreate(outfile,'longitude','Dimensions', {'nlon',size(fluxes,1)},'Datatype','double');
        nccreate(outfile,'latitude','Dimensions', {'nlat',size(fluxes,2)},'Datatype','double');
        ncwrite(outfile,'fluxes',fluxes);
        ncwrite(outfile,'latitude',lat');
        ncwrite(outfile,'longitude',lon);
        ncwriteatt(outfile,'fluxes','unit','micromol m-2 s-1 (averaged over the case study time period)');

	disp('Vector of normed residuals');
	disp(resvec);

end; % End of if statement


%--------------------%
% END OF SECTION     %
%--------------------%


%----------------------------------------------------------------%
% OPTION 2: ESTIMATE UNCERTAINTIES USING REDUCED RANK APPROACH   %
%----------------------------------------------------------------%

        if scriptoption == 2;

        disp('Estimating uncertainties using reduced rank approach');

        uncert_est(D, E, X, theta, n, eigmax, H, outpath, selu);

        end; % End of scriptoption if statement


%-----------------------------------------------------------------------------------------------------------------------%
% END OF SCRIPT														%
%-----------------------------------------------------------------------------------------------------------------------%
