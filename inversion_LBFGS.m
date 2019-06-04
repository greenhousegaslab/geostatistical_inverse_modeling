%-----------------------------------------------------------------------------------------------------------------------%
% SCRIPT:  inversion_run												%
% PURPOSE: Run a geostatistical inverse model (GIM), implemented with Lagrange multipliers 				%
% S. Miller, Jan. 8, 2016												%
%															%
%-----------------------------------------------------------------------------------------------------------------------%

%---------------------------%
% NOTES:                    %
%---------------------------%

	% L-BFGS algorithm
	% Here is the link to the L-BFGS algorithm that I've implemented in this code: 
	% https://www.mathworks.com/matlabcentral/fileexchange/23245-fminlbfgs-fast-limited-memory-optimizer
	% Note: there are several codes available that implement L-BFGS in Matlab. 

	% Before running this script, make sure to download the L-BFGS Matlab scripts from the link above.


%----------------------------%
% Initial setup parameters   %
%----------------------------%

	%---------------------------------------------------------%
	% Set the number of iterations for the L-BFGS algorithm   %
	%---------------------------------------------------------%

	maxit = 50;
	
	%-----------------------------------------%
	% **Set the name of the output folder**   %
	%-----------------------------------------%

	% This script will save the estimated fluxes into the following folder:
	outpath = '';
	
	
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
        %! CHANGE THE PARAMETERS BELOW FOR YOUR PARTICULAR SETUP
        theta = [ 16 0.128 189.7 105 ];

        % Display the covariance matrix parameters on screen
        disp('Covariance matrix parameters');
        disp(theta);
	

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
        % The footprint matrix should have m columns (where m are the total number of fluxes to be estimated).
        % m = q x r
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

        % The X matrix contains spatial and temporal patterns that will be used in the deterministic model of the inversion.
        % The X matrix should have m rows. These rows should correspond to the columns of H.
        % The X matrix will have p columns. Each column contains land surface, environmental, or emissions information that we want to use in the deterministic model.

        % A Bayesian inversion or GIM always assumes that the prior or deterministic model is unbiased relative to the observations. If the prior model is biased relative
        % to the observations, then we have violated the statistical assumptions of the inversion.
        % The GIM has a built-in mechanism to ensure that the deterministic model is not biased relative to the observations.
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
        Dinv = inv(D);


%------------------------%
% Create the R matrix    %
%------------------------%

        disp('Create R');

        % No need to edit this section, unless you want a more complicated setup for R.
        R = (theta(1).^2) .* ones(n,1);
        R = spdiags(R,0,n,n);

        % Define the sizes of the different matrices
        p   = size(X,2);
        m1  = size(E,1);
        ntimes = size(D,1);
        m   = ntimes .* m1;


%---------------------------------------------------%
% Pass the columns of X through the forward model   %
%---------------------------------------------------%

        %! Note: Edit this section to match the actual format of H in your problem.

        disp('Calculate HX');
        HX = zeros(n,p);
        for j = 1:ntimes;
        load(strcat(Hpath,'H_',num2str(j),'.mat'));
        sel = (m1.*(j-1)+1):(j.*m1);
        HX = HX + H*X(sel,:);
        clear H;
        end;


%------------------------------------------------------%
% Create the initial guess for the L-BFGS-B algorithm  %
%------------------------------------------------------%

        disp('Create an initial guess for the L-BFGS-B algorithm');

        % No need to edit this section

        % The L-BFGS algorithm requires an initial guess of the fluxes/emissions. It will then iterate toward the solution.
        % The initial guess can be very important, and a more accurate initial guess is always better than a poor initial guess.
        % Here, I make a guess for beta and use X*beta as the initial guess.
	% Note: In the article Miller et al. (GMD), we use an initial guess of zero. In that paper, we were more interested
	% in understanding the behavior of different algorithms than optimizing the speed of convergence.

        beta1 = (HX' * HX) \ (HX' * Z);
        shat0 = X * beta1;

	%----------------------------------------%
	% Take the symmetric sq root of D and E  %
	%----------------------------------------%

	CD = sqrtm(D);
	CE = sqrtm(E);
	invCD = inv(CD);
	invCE = inv(CE);

	%-------------------------------------------%
	% Transform shat0: shat0* = Q^(-1/2)*shat0  %
	%-------------------------------------------%

	stemp = [];
 
 	ntimes = size(D,1);
         for j = 1:ntimes;
         Qx1   = zeros(m1,1);
                 for i = 1:ntimes;
                 sel = (m1.*(i-1)+1):(i.*m1);
                Qx1 = Qx1 + shat0(sel) .* invCD(j,i);
                end; % End of i loop
        temp =  invCE * Qx1;
        stemp   =  [stemp; temp];
        end; % End of i loop
        clear Qx1 temp;

	shat0 = stemp;


%-------------------------------------------------%
% Pre-calculate matrix products where possible    %
%-------------------------------------------------%

	% Calculate inv(Q) * X
	% We'll refer to this matrix product as varible "B" from now on
	% This step is really slow and is the time-limiting step in the inverse model

	disp('Pre-calculate matrix products where possible');
	
	% B = inv(Q) * X
	B = [];
	p = size(X,2);
	m1  = size(E,1);
	ntimes = size(D,1);
	m = ntimes .* m1;

	for j = 1:ntimes;
	B1 = zeros(m1,size(X,2));
		for i = 1:size(Dinv,1);
		sel = (m1.*(i-1)+1):(i.*m1);
		B1 = B1 + X(sel,:) .* Dinv(j,i);
		end; % End of i loop
	temp = Einv * B1;
	B = [B; temp];
	end; % End of j loop
	clear B1 temp;

	% To save time, one can save out the object B and read it in for future 
	% inverse modeling simulations.
	% save(strcat(outpath,'B.mat'),'B');

	
%-----------------------------------%
% Estimate the fluxes using L-BFGS  %
%-----------------------------------%

        % Set the time counter
        disp('Time at the initiation of the L-BFGS algorithm:');
        disp(clock);


%-------------------------------------------------------------%
% Create function handles for the cost function and gradient  %
%-------------------------------------------------------------%


        %-------------------------------%
        % Calculate inv(sqrtm(Q)) * X   %
        %-------------------------------%

        A = [];

        for j = 1:ntimes;
        A1 = zeros(m1,size(X,2));
                for i = 1:ntimes;
                sel = (m1.*(i-1)+1):(i.*m1);
                A1 = A1 + X(sel,:) .* invCD(j,i);
                end; % End of i loop
        temp = invCE * A1;
        A = [A; temp];
        end; % End of i loop
        clear A1 temp;

	f1 = @(shat) cost_gradient_fun_transform(Z, R, X, A, B, Dinv, Einv, CD, CE, Hpath, shat);

        % Create an empty flux estimate
        shat = [];


%-----------------------------------------%
% Set options for the L-BFGS-B algorithm  %
%-----------------------------------------%

        options = struct('HessUpdate','lbfgs','GradObj','on','Display','iter','MaxIter',maxit,'GradConstr',false);


%-----------------------%
% Run the algorithm     %
%-----------------------%

        disp('Run the L-BFGS-B algorithm');

        [shat,costfun,exitflag,gradient] = fminlbfgs(f1,shat0,options);


%--------------------------------------------%
% Transform fluxes back to normal space      %
%--------------------------------------------%

	stemp = [];
 
 	ntimes = size(D,1);
         for j = 1:ntimes;
         Qx1   = zeros(m1,1);
                 for i = 1:ntimes;
                 sel = (m1.*(i-1)+1):(i.*m1);
                Qx1 = Qx1 + shat(sel) .* CD(j,i);
                end; % End of i loop
        temp =  CE * Qx1;
        stemp   =  [stemp; temp];
        end; % End of i loop
        clear Qx1 temp;

	shat = stemp;

        % Print out time information
        disp('Time at the end of the L-BFGS algorithm:');
        disp(clock);
	
	
%------------------------------%
% Write the outputs to file    %
%------------------------------%

        disp('Writing outputs to file');
        dlmwrite(outname,full(shat),',');
	save(strcat(outpath,'fluxes_LBFGS',num2str(maxit),'.mat'),'shat');
	
        disp('Outputs written to file');
        disp(outname);	


%-----------------------------------------------------------------------------------------------------------------------%
% END OF SCRIPT														%
%-----------------------------------------------------------------------------------------------------------------------%
