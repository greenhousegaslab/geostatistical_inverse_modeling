%---------------------------------------------------------------------------------------------------------------%
% FUNCTION: model_selection_run.m										%
% PURPOSE: Run a BIC-based model selection.									%
% S. Miller, Aug. 9, 2016											%
%	Updated Aug 28, 2025											%
%														%
%---------------------------------------------------------------------------------------------------------------%

	%-----------%
	% NOTES:    %
	%-----------%

	% I've designed this script to output the model that has the lowest BIC score and also has non-negative coefficients.
	% 	In the case of CH4, negative coefficients are usually physically unrealistic, and I throw out models that produce
	% 	negative coefficients. You may want to change this setup, however, based upon the problem of interest. For example,
	% 	perhaps negative coefficients would be unrealistic for some auxiliary variables but not for others.
	%	E.g., CO2 fluxes from the biosphere could have a negative relationship with some environmental variables.
	%	(I.e., higher temperatures yield more negative CO2 fluxes.)

	% Any heading that begins with "**" requires a user input.


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

	% Each column of the X matrix should be a candidate auxiliary variable.

        X = [ READ IN THE X MATRIX HERE ];

	% Record the dimensions of X
	p = size(X,2);
	m = size(X,1);

	% Add an intercept term to the X matrix
	% The intercept term is a column of ones.
	X = [ ones(m,1) X ];	


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

        sigmaQ = theta(2) .* ones(m,1);

        D = (sigmaQ*sigmaQ') .* D;


%------------------------%
% Create the R matrix    %
%------------------------%

        disp('Create R');

        % No need to edit this section, unless you want a more complicated setup for R.
        R = (theta(1).^2) .* ones(n,1);
        R = spdiags(R,0,n,n);


%---------------------------------------------------%
% Pass the columns of X through the forward model   %
%---------------------------------------------------%

	% Call the "metvecH" function to create an H matrix object.
        ntimes = size(D,1);

	H = matvecH(ntimes,Hpath);
	HX = H*X;


%---------------------------%
% Calculate psi and psi^-1  %
%---------------------------%

	disp('Compute psi and the inverse of psi');

        % Compute the psi matrix
        psi = HQHR(R, D, E, H);

	% Calculate the determinant of psi
        detpsi = 2.*sum(log(diag(chol(psi))));

        % Take the inverse of psi
        %!! This step might be really really slow if n is large. 
        % If so, contact me, and we can chat about more efficient (but more involved) computational approaches.
        psiinv = psi \ speye(n,n);
        clear psi;


%----------------------------------%
% Create the A, B, and C matrices  %
%----------------------------------%

        % For more detail on these matrices, see Vineet Yadav's paper on model selection.

        A = Z' * psiinv * Z;
        B = Z' * (psiinv * HX);
        C = HX' * psiinv * HX;


%----------------------------%
% Check for colinearity      %
%----------------------------%

	% The intercept term is excluded from the colinearity calculations here.

	%-------%
	% VIF   %
	%-------%

        p0 = size(C,1) - p;
        VIF = abs(diag(inv(corrcoef( HX(:,(p0+1):(p0+p)) ))));

        disp('Variance inflation factors');
	disp('(A VIF greater than 10 incidates colinearity)');
        disp(VIF);

        clear HX;

	%---------------------%
	% Condition number:   %
	%---------------------%

	condnum = cond(C);
	
	disp('Condition number');
	disp('(A condition number greater than 30 can indicate colinearity)');
	disp(num2str(condnum));	


%----------------------------------------------------------%
% Calculate a BIC score for each possible candidate model  %
%----------------------------------------------------------%

	disp('Calculate BIC scores for all possible candidate models');

	[ selmodel beta1 ] = bic(Z, A, B, C, p, detpsi);

	% Display the outputs on screen
	disp('Column index of the selected best model');
	disp(num2str(selmodel));
	

%---------------------------------------------------------------------------------------------------------------%
% END OF SCRIPT													%
%---------------------------------------------------------------------------------------------------------------%
