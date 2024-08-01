%-----------------------------------------------------------------------------------------------%
% FUNCTION: uncert_est.m									%
% PURPOSE:  Approximate the posterior uncertainties using a reduced rank approach.		%
% S. Miller, Nov. 16, 2018									%
%	Updated July 24, 2024									%
%												%
%-----------------------------------------------------------------------------------------------%


%-----------%
% NOTES:    %
%-----------%

	% Note: The cholesky decomposition in Matlab is faster than the matrix square root function in Matlab,
	% so I use the former here.

	% If you don't use enough eigvenvalues/vectors, then inv(V3) can produce negative values.
	% This will result in an imaginary 95% confidence interval.

	% Note: In this script, we assume that the fluxes are being estimated at a 3-hourly time resolution 
	% and have units of micromol m-2 s-1. Furthermore, we assume that the results are being reported 
	% in units of petagrams carbon. If that is not the case, make sure to edit the units conversion 
	% where noted in the script below.

	% Furthermore, the most computationally intensive part of the script is the estimation of eigenvalues
	% and vectors. The script, as written below, will write the eigenvalues and vectors to a .mat file.
	% Every time the script is run, it will check to see if a file with eigenvalues and vectors exists.
	% If so, it will read in that file instead of calculating new eigenvalues and vectors.


%-------------------%
% BEGIN FUNCTION    %
%-------------------%

	function [ ] = uncert_est(D, E, X, theta, n, eigmax, H, outpath, selu);

	% FUNCTION INPUTS:
	% D,E: 		Components of the Q matrix
	% X:		Matrix that contains the auxiliary variables.
	% theta:	Vector with parameters that define the covariance matrices.
	% n:		Number of observations in the inverse problem.
	% eigmax:	Maximum number of eigenvalues and vectors to calculate.
	% H:		Object of class matvecH that points to the H matrix strips.
	% outpath:	Location where the outputs should be saved.
	% selu:		Vector that contains the area of the model grid boxes that should
        %		be included in the uncertainty calculations. Fill the vector with 
	%		zeros for grid boxes that should be excluded from the uncertainty calculations. 

	% FUNCTION OUTPUTS:
	% None (results are printed on screen).


%------------------------------------%
% Optional user-defined parameters   %
%------------------------------------%

	% File containing the eigenvalues/vectors (for option 1)
	% File where the eigenvalues/vectors should be written (for option 2)
	eigfile = strcat(outpath,'eigmat_',num2str(eigmax),'.mat');

	% This script will search for the file containing the eigenvalues. It'll use that file if it exists.
	% If not, it will create that file.


%-------------------------------%
% Essential script parameters   %
%-------------------------------%

        ntimes = size(D,1);
        m = size(X,1);
        m1 = m ./ ntimes; % Divide m by the total number of time periods in the inversion


%----------------------------------------------%
% Find the Cholesky decomposition of D and E   %
%----------------------------------------------%

	CD   = chol(D)'; % Transpose has better numerical stability in Matlab
	CE   = chol(E)';
	Dinv = inv(D);
	Einv = inv(E);

	Qchol = kronMat(CD,CE);
	Qinv  = kronMat(Dine,Einv);


%---------------------------------------------%
% Check whether the eigenvalue file exists    %
%---------------------------------------------%

	eigoption = exist(eigfile);

%-----------------------------------------------%
% OPTION 1: Compute the eigenvalues/vectors     %
%-----------------------------------------------%

	if eigoption == 0;

	% Make a counter to tally the number of times "eigs" calls "eigfun"
	counter = 0;
	save(strcat(Hpath,'counter_',num2str(eigmax),'.mat'),'counter');

	f1 = @(x1) eigfun(Qchol, theta, n, H, x1);

	disp('Number of eigenvalues/vectors to calculate');
	disp(num2str(eigmax));


	%-----------------------------------------------%
	% Option (a): Run the eigs function in Matlab   %
	%-----------------------------------------------%
	
	% options = struct('IsFunctionSymmetric',1);
	% [eigvec,eigval] = eigs(f1,m,eigmax,'LM',options);

	%---------------------------------------------%
	% Option (b): Run randomized eigen estimator  %
	%---------------------------------------------%

	[eigvec,eigval] = randeigdecomp(Qchol, theta, n, H, eigmax);

	% Save the calculated eigenvalues and vectors to file
	save(eigfile,'eigvec','eigval','-v7.3');

	end; % End of if statement


%------------------------------------------------%
% OPTION 2: Read eigenvalues/vectors from file   %
%------------------------------------------------%

	if eigoption == 2;

	load(eigfile);

	end; % End of eigoption if statement

	disp('Eigenvalues');
	disp(num2str(diag(eigval)));

	eigval = eigval./(eigval+1);


%----------------------%
% Compute inv(Q)*X     %
%----------------------%

        % Calculate inv(Q) * X

        % B = inv(Q) * X
	B = Qinv * X;


%-----------------------------------------------------------------%
% Calculate the uncertainty in the total flux from a given region %
%-----------------------------------------------------------------%

	% Repeat the unit conversion vector based upon the number of time periods in the inverse model
	agvec = repmat(selu,ntimes,1);


%----------------------------------------%
% Calculate the posterior uncertainties  %
%----------------------------------------%

	% V_shat = V1 + V2*V3*V2
	% V1 = Q - Q^(1/2)*B*A*B'*Q^(1/2)
	% V2 = V1*inv(Q)*X
	% V3 = inv(X'*inv(Q)*X - [inv(Q)*X]'*V1*[inv(Q)*X]) 


        %----------------------------------%
        % Mutiply V1*(aggregation vector)  %
        %----------------------------------%

        % Multiply Q by aggregation vector
        Qag = Q * agvec;

        % Finish calculating (agvec)' * Q * agvec
        V1a = agvec' * Qag;

	disp('V1');
	disp(num2str(V1a));


	%-------------------------%
	% Finish calculating V1   %
	%-------------------------%

        % Multiply QBABQ by the aggregation vector
        [ QBABQx ] = QBABQ(CD, CE, eigvec, eigval, agvec);
        V1b = agvec' * QBABQx;
	clear QBABQx;


	%---------------%
	% Compute V2    %
	%---------------%

	% V2 = V1*inv(Q)*X
	% V2 = Q*inv(Q)*X - Q^(1/2)*B*A*B'*Q^(1/2)*inv(Q)*X
	% V2 = X - Q^(1/2)*B*A*B'*Q^(1/2)*inv(Q)*X (Could simplify Q, but unlikely to save compute time)

	[ QBABQx ] = QBABQ(Qchol, eigvec, eigval, QX);

	% Subtract from X to get V2
	V2 = X - QBABQx;	
	clear QBABQx;


	%----------------%
	% Compute V3     %
	%----------------%

	% V3 = inv(X'*inv(Q)*X - [inv(Q)*X]'*V1*[inv(Q)*X]) 
	% V3 = inv(X'*inv(Q)*X - [inv(Q)*X]' * V2)

	V3 = X' * QX - QX' * V2; 
	V3 = inv(V3);


	%------------------------------------------%
	% Multiply V2*V3*V2*(aggregation vector)   %
	%------------------------------------------%

	V2V3V2 = agvec' * (V2 * (V3 * (V2' * agvec)));


	%------------------------------------------%
	% Compute the final uncertainty estimate   %
	%------------------------------------------%

	uncert = V1a - V1b + V2V3V2;
	uncert = sqrt(uncert);

	% Note: the lines below will convert the fluxes from units of micromol / m-2 
	% per 3-hourly time period of the inverse model to units of petagrams of carbon
	% summed over the entire time period of the inverse model.
	% Edit the code below if you are using different units in the inverse model.
	% Convert from micromol to Pg C
        % Convert micromol to mol
        uncert = uncert ./ 1e6;

        % Convert mol to gram C
	uncert = uncert .* 12;

        % Convert gram to petagram
	uncert = uncert .* 1e-15;

	% Convert to 95% confidence interval
	% Uncertainties are in Pg C
	uncert = 1.96 .* uncert;

	disp('95% confidence interval equals the best')
	disp('estimate plus or minus this number:');
	disp(num2str(uncert));


%-----------------------------------------------------------------------------------------------%
% END OF FUNCTION										%				
%-----------------------------------------------------------------------------------------------%

