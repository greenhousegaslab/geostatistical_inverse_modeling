%-----------------------------------------------------------------------------------------------------------------------%
% SCRIPT:  inversion_dual_eq												%
% PURPOSE: Launch a geostatistical inverse model (GIM).									%
%	   This script will use the analytical equations in Michalak et al. (2004) to estimate the fluxes.		%
%	   I.e., this script does not estimate the fluxes iteratively.							%
% S. Miller, Jan. 8, 2016												%
%															%
%-----------------------------------------------------------------------------------------------------------------------%

%---------------------------%
% Summary of this script    %
%---------------------------%

	% In this script, I implement an inverse model using the minimum residual method to iterate toward the solution.
	% I can then compare the accuracy of the flux estimate from this script against the true analytical solution 

	% Note that the preconditioner didn't really help for the tests that I did. I tried both 100 total eigenvectors and 400 total eigenvectors.
	% The convergence results were pretty similar for both.


%---------------------------------%
% REQUIRED USER-DEFINED OPTIONS   %
%---------------------------------%

	% Set what you want to calculate (e.g., best estimate, uncertainties). See description of options below.
	scriptoption = 4;

	% Set path where outputs should be saved
	outpath = '/home-4/smill191@jhu.edu/work/smiller/GIM_test/red_rank/';


%--------------------------------------------------------------------%
% Options 1 and 2: estimate fluxes with or without a preconditioner  %
%--------------------------------------------------------------------%

        if scriptoption < 3;

        disp('Calculate the posterior best estimate of the fluxes');

        % Set the number of iterations for the minimum residual algorithm
        % Read in this information from the unix enviroment
        % (Another option is to define the variable explicitely here)
        maxit = str2num(getenv('maxit'));
        disp('maxit');
        disp(num2str(maxit));

        % This script will save the estimated fluxes into the following file (csv file):
        outfile = strcat(outpath,'fluxes_minresid',num2str(maxit),'.csv');
        weightfile = strcat(outpath,'weights_',num2str(maxit),'.mat');

	% Read in the analytical solution, and I can compare the accuracy of the estimate as it converges
	shata = dlmread('/home-4/smill191@jhu.edu/work/smiller/GIM_test/fluxes.csv');

	% Delete any existing information on algorithm convergence
	% if maxit == 1;
        % unix(strcat(outpath,'convergence_minresid.mat'));
        % convergence = [];
        % save(strcat(outpath,'convergence_minresid.mat'),'convergence');
	% end;	

	end; % End of scriptoption if statement


%-------------------------------------------------------------------%
% Option 3: estimate uncertainties using conditional realizations   %
%-------------------------------------------------------------------%

        if scriptoption == 3;

	disp('Generate conditional realizations');

	% Set the maximum number of iterations for the minres algorithm
	maxit = 50;

	% Set the maximum number of conditional realizations
	maxcond = 500;

	end; % End of scriptoption if statement


%--------------------------------------------------------------%
% Option 4: estimate uncertainties with reduced rank approach  %
%--------------------------------------------------------------%

	if scriptoption == 4;

        disp('Estimating uncertainties using reduced rank approach');

        eigmax = str2num(getenv('eigmax'));
        disp('eigmax');
        disp(num2str(eigmax));
	
	end; 


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


%---------------------------------%
% **Read in the observations**    %
%---------------------------------%

	disp('Read in the observation vector');

	% Read in the vector of observations here.
	% Subtract off the boundary condition before reading in the observations here.  
	% We'll refer to this vector as "Z".
	% This vector has length (n x 1) where n are the number of observations

        % load('/home-4/smill191@jhu.edu/work/smiller/GIM_test/z.mat','z')
	% z = sum(z,2);
        % Z = z + normrnd(0,theta(1),[length(z),1]);
	% save('/home-4/smill191@jhu.edu/work/smiller/GIM_test/Z.mat','Z');
	load('/home-4/smill191@jhu.edu/work/smiller/GIM_test/Z.mat');


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
	
	% NOTE: I've re-written the inverse modeling code so that we can ignore this section. The code below will read in smaller pieces of H
	% and not read in the entire H matrix at once.

	Hpath = '~/work/smiller/data/OCO2_H_matrix/';


%----------------------------%
% **Create the X matrix**    %
%----------------------------%

	disp('Create the X matrix');

	% The X matrix contains spatial and temporal patterns that will be used in the prior model of the inversion.
	% In a GIM, we usually refer to this prior model as the "deterministic model"
	% The X matrix should have (q x r) rows. These rows should correspond to the columns of H.
	% The X matrix will have p columns. Each column contains land surface, environmental, or emissions information that we want to use in the deterministic model.
	
	% Options for X:
	% OPTION 1: In the simplest case, X could be a vector of ones. In this case, X will have p=1 columns. This setup is similar to what Kim Mueller had done in her 2008 CO2 inversion paper. 
	%	If you set up the X matrix in this way, you are making as few assumptions about the fluxes as possible. Any spatial/temporal patterns in the fluxes will come directly
	%	from the observations and not from an emissions inventory or bottom-up emissions model.
	% OPTION 2: One could include an existing emissions inventory in the X matrix (e.g., the EDGAR inventory). In this setup, X would have p=2 columns. The first column would be a 
	%	vector of ones, and the 2nd column would be the chosen emissions inventory. 
	%	In a GIM, it is customary to include a vector of ones in X, irrespective of 
	%	what else gets included in X. This vector of ones is analagous to the intercept term in a regression. In a GIM, this vector of ones will account for any patterns in the fluxes
	%	that aren't captured by other columns of the X matrix. 
	% OPTION 3: X could have a more complex format with land surface data, meteorological reanalysis, etc. Sharon Gourdji's 2008 and 2012 papers provide good examples. 

	tic;	
	n = length(Z);

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
	
	disp(toc);


%-------------------------%
% Create the E matrix     %
%-------------------------%

	disp('Create the E and D matrices');

	tic;

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

	% load('~/data/smiller/outputs/xiaoling/distmat.mat');	
	load('/home-4/smill191@jhu.edu/work/smiller/GIM_test/distmat.mat');	


	%------------%
	% Create E   %
	%------------%

	% No need to edit this section
	% Note: I use a spherical covariance model because of its sparse properties.

	% Multiply the distance matrix by the covariance model
	% Exponential covariance model
	% E = exp(-1 .* deltamat ./ theta(3));

	% Spherical covariance model
	E                      = 1 - 1.5 .* (deltamat ./theta(3))  + 0.5 .* (deltamat.^3 ./ theta(3).^3);
	E(deltamat > theta(3)) = 0;

	Einv = inv(E);
	
	
%----------------------------%
% Create the D matrix        %
%----------------------------%

	% The D matrix describes the time between each time period in the inversion
	% It will have dimensions (q x q)
	% The D matrix is usually easier to set up than the E matrix, and I have included sample code below
	
	%-----------------------------------%
	% **Create time distance matrix**   %
	%-----------------------------------%

    	days = [];
    	for i = 1:328
        for j = 1:328
            if rem(abs(i-j),8) == 0
                days(i,j) = abs(i-j)/8;
            else
                days(i,j) = 10^6;
            end
        end
    	end

	%------------%
	% Create D   %
	%------------%

	% No need to edit this section
	% Exponential model
	% D      = exp(-1 .* days ./ theta(4));
	
	% Spherical covariance model
	D = 1 - 1.5 .* (days   ./theta(4))  + 0.5 .* (days.^3   ./ theta(4).^3);
	D(days   > theta(4)) = 0;

	disp(toc);


%----------------------------%
% Create the sigmaQ vector   %
%----------------------------%

	% This section defines the diagonal elements of the Q matrix.
	% There is no need to edit this section if you want Q to have constant diagonal elements.
	% If you want Q to have varying diagonal elements, then you will need to edit this section.

	m = size(X,1);
	m1= size(D,1);	

	sigmaQ = theta(2) .* ones(m1,1);

        D = (sigmaQ*sigmaQ') .* D;
        Dinv = inv(D);


%------------------------%
% Create the R matrix    %
%------------------------%
		
	disp('Create R');
	tic; 
	
	% No need to edit this section, unless you want a more complicated setup for R.
	R = (theta(1).^2) .* ones(n,1);
	R = spdiags(R,0,n,n);

	disp(toc);


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
        tic;
        HX = zeros(n,p);
        for j = 1:size(D,1);
        load(strcat(Hpath,'H_',num2str(j),'.mat'));
        sel = (m1.*(j-1)+1):(j.*m1);
        HX = HX + H*X(sel,:);
        clear H;
        end;
        disp(toc);


%-------------------------------------------------------------------------%
% OPTION 1: USE A PRECONDITIONER WITH MIN RES TO CALCULATE BEST ESTIMATE  %
%-------------------------------------------------------------------------%

	%------------------------------------------%
	% Check whether the fluxes already exist   %
	%------------------------------------------%

	fflag = 0;

	if scriptoption < 3;

	fflag = exist(outfile) + exist(weightfile);

	if fflag == 4;
	disp('Existing flux estimate found. Read in this estimate.');
	shat = csvread(outfile);
	load(weightfile);
	end;

	end;


        % Details on the Matlab minimum residual function:
        % x = minres(A,b) attempts to find a minimum norm residual solution x to the system of linear equations A*x=b
        % A can also be a function handle that returns A*x

	if scriptoption == 1 & fflag < 4;

	disp('Pre-calculate items for the preconditioner');


	%--------------------------------------------------%
	% Check whether the preconditioner already exists  % 
	%--------------------------------------------------%

	prefile = strcat(outpath,"preconditioner.mat");

	if exist(prefile)==2;

	disp('Preconditioner already exists. Read from file');
	load(prefile);

	else;

	disp('Preconditioner does not exist. Calculate.');

        %------------------------------------------------------%
        % Pre-calculate items required by the preconditioner   %
        %------------------------------------------------------%

	%----------------------------------------%
	% Compute the eigenvalues/vectors of Q   %
	%----------------------------------------%

        % The eigenvectors/values of Q should be the Kronecker product of the eigenvectors of D and E
        % and the Kronecker product of the eigenvalues of D and E

	% Number of eigenvectors to generate from D and E
	% The greater the number of eigenvectors, the better (but more computationally intensive) the preconditioner
	neigs = 20;

        % Compute the first 10 eigenvectors/values of D
        % Use the "eigs" function in Matlab
        [VD,DD] = eigs(D,neigs);

        % Compute the first 10 eigenvectors/values of E
        [VE,DE] = eigs(E,neigs);

        % Use the Kronecker product to generate full eigvenvectors
        Vall = kron(VD,VE);

        % Use the Kronecker prodcuct to generate the full eigenvalues
        Dall = kron(DD,DE);


	%----------------------------%
	% Multiply eigvectors by H   %
	%----------------------------%

	% Step 2 in Saibaba et al. 2012

        n = size(R,1);
        m1  = size(E,1);

        % disp('Calculate HV');
        HVall = zeros(n,size(Vall,2));
        for j = 1:size(D,1);
        load(strcat(Hpath,'H_',num2str(j),'.mat'));
        sel = (m1.*(j-1)+1):(j.*m1);
        HVall = HVall + H*Vall(sel,:);
        clear H;
        end;

	%---------------%
	% Calculate M   %
	%---------------%

	% Step 2 in Saibaba et al. 2012
        % This matrix will have dimensions n x 100

        sqR = (theta(1).^2).^(-0.5) .* ones(n,1);
        sqR = spdiags(sqR,0,n,n);
        % sqR = R.^(-0.5);

        M = sqR * HVall * (Dall .^ 0.5);

	%-------------------------------------%
	% Singular value decomposition of M   %
	%-------------------------------------%

	% Step 3 in Saibaba et al 2012

        % Note: This step takes a lot of memory
	% The econ option in svd will produce a square Sigma1 matrix
        [Um,Sigma1, Vm] = svd(M,'econ');
        clear Vm M;

        % I'm going to chop down the dimensions of Um.
        % I'm not sure if this step is correct, but it's the best I know how to do.
        % Um = Um(:,1:(neigs.^2));


	%----------------------------------------------%
	% Calculate D using Sherman Morrison update    %
	%----------------------------------------------%

	% Step 4 in Saibaba et al. 2012

        Sigma1 = diag(Sigma1);
        Dr     = diag(Sigma1.^2 ./ (1 + Sigma1.^2));

	%--------------------------------%
	% Write preconditioner to file   %
	%--------------------------------%

	save(prefile,'Um','Dr','sqR');

	end; % End of if else statement

	%-------------------------------------%
	% Launch the minimum residual method  %
	%-------------------------------------%

	% Set the function options
	tol   = []; % Tolerance
	% maxit = 3; % Max number of iterations
	% M     = []; % A pre-conditioner for the system

	% Set a first guess
	% Note: if no first guess is specified, then the first guess is automatically set at zero.	
	
	% Create a function handle
	f1 = @(weights) Ax(R, X, HX, D, E, Hpath, weights);
	f2 = @(x1) preconditioner(R, sqR, Um, Dr, HX, x1);;

	% Minimum residual method with preconditioner	
	[weights1, flag,relres,iter,resvec] = minres(f1,b,tol,maxit,f2);
	% [weights1] = minres(f1,b,tol,maxit,f2);

	disp('Residual vector');
	disp(num2str(resvec));

	end; % End of if statement

%-----------------------------------%
% END OF OPTION 1 		    %
%-----------------------------------%

%-------------------------------------------------------------%
% OPTION 2: ESTIMATE FLUXES USING MIN RES, NO PRECONDITIONER  %
%-------------------------------------------------------------%

        if scriptoption == 2 & fflag < 4;

	% Minimum residual method without preconditioner
	disp('Estimate weights using minimum residual method');
	tol = []; % Tolerance
	f1 = @(weights) Ax(R, X, HX, D, E, Hpath, weights);
        [weights1, flag,relres,iter,resvec] = minres(f1,b,tol,maxit);

	% Conjugate gradients square method
        % [weights1, flag,relres,iter,resvec] = cgs(f1,b,tol,maxit,M);

	% Biconjugate gradients method
        % Note: this function requires a transpose option in Ax (see Matlab web page for more details)
	% [weights1, flag,relres,iter,resvec] = bicg(f1,b,tol,maxit,M);

	end; % End of option 2

%-----------------------------------%
% END OF OPTION 1                   %
%-----------------------------------%

%---------------------------------------------------------------------%
% OPTIONS 1 AND 2: COMPUTE FLUXES FROM MINRES OUTPUT  		      %
%---------------------------------------------------------------------%

	if scriptoption < 3; 

	if fflag < 4;

	shat = weights_to_fluxes(weights1,Hpath,D,E,X,p);


	%------------------------------%
	% Write the outputs to file    %
	%------------------------------%

        disp('Writing outputs to file');
        dlmwrite(outfile,full(shat),',');
	save(weightfile,'weights1');
	
        disp('Outputs written to file');
        disp(outfile);	

	disp('Vector of normed residuals');
	disp(resvec);

	end; % End of fflag if statement


	%----------------------------------------------------------------------------%
	% Calculate the root mean squared error relative to the analytical solution  %
	%----------------------------------------------------------------------------%

        % Calculate the root mean squared error
        rmse = sqrt(mean((shat - shata).^2));

        disp('Current root mean squared error');
        disp('(relative to analytical solution)');
        disp(num2str(rmse));


	%-------------------------------------------------------------%
	% Calculate the value of the inverse modeling cost function   %
	%-------------------------------------------------------------%

        % Calculate Hs
        Hs = zeros(n,1);
        for j = 1:size(D,1);
        load(strcat(Hpath,'H_',num2str(j),'.mat'));
        sel   = (m1.*(j-1)+1):(j.*m1);
        Hs = Hs + H*shat(sel,:);
        clear H;
        end;

        % Calculate Q^-1 (s - Xbeta)
        betas = weights1((size(weights1,1)-p+1):size(weights1,1));;
	sXbeta  = (shat - X*betas);
        QsXbeta = [];
        for j = 1:size(Dinv,1);
        A1 = zeros(m1,1);
                for i = 1:size(Dinv,1);
                sel = (m1.*(i-1)+1):(i.*m1);
                A1 = A1 + sXbeta(sel) .* Dinv(j,i);
                end; % End of i loop
        temp = Einv * A1;
        QsXbeta = [QsXbeta; temp];
        end; % End of i loop
        clear A1 temp;

        % Calculate the cost function
        costfun1 = (Z - Hs)' * (R \ (Z - Hs));
        costfun2 = sXbeta' * QsXbeta;
        costfun  = costfun1 + costfun2;

        disp('Cost function value');
        % disp(num2str(costfun1));
        % disp(num2str(costfun2));
        disp(num2str(costfun1 + costfun2));


	%---------------------------------------%
	% Print out the minimum residual norm   %
	%---------------------------------------%

        % norm(b-A*x0)
	A = Ax(R, X, HX, D, E, Hpath, weights1);
        b = [ Z ; zeros(p,1)];
        mrnorm = norm(b - A);

        disp('Minimum residual norm');
        disp(num2str(mrnorm));

        % Save these metrics to a file
        % load(strcat(outpath,'convergence_minresid.mat'));
        % temp = [maxit costfun rmse mrnorm];
	convergence = [maxit costfun rmse];

        % convergence = [convergence; temp];
        save(strcat(outpath,'convergence_minresid_',num2str(maxit),'.mat'),'convergence');

	end; % End of if statement


%--------------------%
% END OF SECTION     %
%--------------------%


%-------------------------------------------------------------------%
% OPTION 3: ESTIMATE UNCERTAINTIES USING CONDITIONAL REALIZATIONS   %
%-------------------------------------------------------------------%

	if scriptoption == 3;

	disp('Begin creating conditional realizations');

	% Create an empty output object
	% Note: I'll only save out the time-averaged flux to save on memory
	condreals = zeros(m1,maxcond);

	% Find the cholesky decomposition of D and E
        CD = chol(D);
        CE = chol(E);

	% Initiate a parallel loop
	parpool(12);
	parfor k = 1:maxcond;
	disp(num2str(k));

	tic;

	% Generate random numbers
	% Create a new random seed
	rng(k);
	% Generate random numbers
	u1 = randn(n,1);
	u2 = randn(m,1);

	% Create random numbers that follow N(0,R)
	u1 = theta(1) .* u1;

	% Create random numbers that follow N(0,Q)
        Qx2 = [];

        for j = 1:ntimes;
        Qx1   = zeros(m1,1);
                for i = 1:ntimes;
                sel = (m1.*(i-1)+1):(i.*m1);
                Qx1 = Qx1 + u2(sel) .* CD(j,i);
                end; % End of i loop
        temp =  CE * Qx1;
        Qx2   =  [Qx2; temp];
        end; % End of i loop
        Qx1 = []; temp = []; 
	u2   = Qx2;

	% Calcualte H*u2
        Hu2 = zeros(n,1);
        for j = 1:size(D,1);
        H = load(strcat(Hpath,'H_',num2str(j),'.mat'));
        H = H.H;
	sel   = (m1.*(j-1)+1):(j.*m1);
        Hu2 = Hu2 + H*u2(sel,:);
        H = [];
        end;
	
	% Generate the conditional realization
	b = [ (Z - Hu2 + u1); zeros(p,1)];
        tol = []; % Tolerance
	f1 = @(weights) Ax(R, X, HX, D, E, Hpath, weights);
       [weights1, flag,relres,iter,resvec] = minres(f1,b,tol,maxit);

	disp('Number of iterations required to converge');
	disp(num2str(iter));

	% Create the flux estimate using the weights
	condreal = weights_to_fluxes(weights1,Hpath,D,E,X,p);
	condreal = condreal + u2;

	% Time-average the flux estimate
	condreal = reshape(condreal,m1,ntimes);
	condreal = mean(condreal,2);

	% Save the outputs
	condreals(:,k) = condreal;
	% condreals = [condreals condreal];

	disp('Time to compute one realization');
	disp(num2str(toc));

	end; % End of parfor loop	


	%---------------------------------------------%
	% Write the conditional realizations to file  %
	%---------------------------------------------%

	disp('Writing conditional realizations to file');
	save(strcat(outpath,'conditional_realizations.mat'),'condreals');
	disp('File written');

	end; % End of scriptoption if statement


%----------------------------------------------------------------%
% OPTION 4: ESTIMATE UNCERTAINTIES USING REDUCED RANK APPROACH   %
%----------------------------------------------------------------%

        if scriptoption == 4;

        disp('Estimating uncertainties using reduced rank approach');

	% Open parallel computing channels
	parpool(5)
	
        uncert_est(D, E, X, theta, n, eigmax, Hpath, outpath);

        end; % End of scriptoption if statement


%-----------------------------------------------------------------------------------------------------------------------%
% END OF SCRIPT														%
%-----------------------------------------------------------------------------------------------------------------------%
