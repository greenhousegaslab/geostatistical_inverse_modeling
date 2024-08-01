%-----------------------------------------------------------------------------------------------------------------------------------------------%
% FUNCTION; cost_gradient_fun															%
% PURPOSE:  Compute the geostatistical inversion cost function. Use the Kronecker product to make these calculations more effeicient.		%
%	     This script will calculate both the cost function and the gradient (when the latter is requested).					%
%	    This script does not use a variable transformation to speed up convergence.								%
%	    This script is useful for inverse problems with bounds but can converge on the solution very slowly.				%
% S. Miller, Sept. 21, 2018															%
%																		%
%-----------------------------------------------------------------------------------------------------------------------------------------------%


	%-------------%
	% NOTES:      %
	%-------------%
		
	
function [ f,g ] = cost_gradient_fun(Z, R, X, B, Qinv, H, shat);

	% FUNCTION INPUTS:
	% Z		Observation vector
	% H		Footprints
	% R		Model-data mismatch matrix
	% X		Design matrix
	% Qinv		Inverse of Q. Formulated as a kronMat object.

	% FUNCTION OUTPUTS:
	% f: 		Value of the cost function for a given guess for the fluxes (shat) 
	% g:		Calculated gradient of the cost function.

	
%----------------------------------------------%
% Set variables that give problem dimensions   %
%----------------------------------------------%

	ntimes = size(Dinv,1);
	n = size(R,1);
	m = size(X,1);
	m1 = m ./ ntimes; % Divide m by the total number of time periods in the inversion


%---------------------------------------%
%***** CALCULATE THE COST FUNCTION *****%
%---------------------------------------%

disp('Calculating the cost function');


%------------------------------------------------------------%
% Pass the current flux estimate through the forward model   %
%------------------------------------------------------------%

	tic;

        %! Note: Edit this section to match the actual format of H in your problem.

        % disp('Calculate Hs');
	Hs = H * shat;


%-------------------------------------------%
% Calculate the likelihood/data component   %
%-------------------------------------------%

	% Fastest formulation
	L1 = sum((Z - Hs).^2./diag(R));

	disp('Cost function value of L1');
	disp(num2str(L1));


%---------------------------------%
% Calculate the prior component   %
%---------------------------------%
	
	% Calculate G * shat
	% G * shat = A - B * inv( X' * B) * B' * shat
	
	%---------------%
	% Calculate A   %
	%---------------%
	
	% A = inv(Q) * shat
	A = Qinv * shat;

	
	%------------------------------%
	% Calculate shat' * G * shat   %
	%------------------------------%
	
	Gshat = A - B * ((X' * B) \ (B' * shat));
	
	L2 = shat' * Gshat;
	clear Gshat;

	disp('Cost function value of L2');
	disp(num2str(L2));
	

%-------------------------------------------------%
% Add up the two components of the cost function  %
%-------------------------------------------------%

	f = L1 + L2;

	disp('Time elapsed for cost function calculations');
	disp(toc);

	disp('Cost function value');
	disp(f);


%---------------------------------------%
%***** CALCULATE THE GRADIENT      *****%
%---------------------------------------%

if ( nargout > 1 );

disp('Calculating the gradient');


%------------------------------------------------------%
% Calculate the observation component of the gradient  %
%------------------------------------------------------%

	% In this section of the code, we want to do the following matrix/vector calculations:
	% L1 = H' * (R \ (Z - Hs));
	% To do these calculations, we'll need to use the GEOS-Chem adjoint

	tic;
	temp = (R \ (Z - Hs));

        %! Note: Edit this section to match the actual format of H in your problem.

	% Calculate  L1 = H' * (R \ (Z - Hs))] 
        % disp('Calculate H * (R \ (Z - Hs))');
	L1 = H' * temp;


%------------------------------------------------%
% Calculate the prior component of the gradient  %
%------------------------------------------------%

	% Calculate G * shat
	% G * shat = A - B * inv( X' * B) * B' * shat
	
	%---------------%
	% Calculate A   %
	%---------------%
	
	% A = inv(Q) * shat
	A = Qinv * shat;
	
	
	%------------------------------%
	% Calculate shat' * G * shat   %
	%------------------------------%
	
	L2 = A - B * ((X' * B) \ (B' * shat));
	

%--------------------------------------------%
% Add up the two components of the gradient  %
%--------------------------------------------%

	g = L2 - L1;

	disp('Time elapsed for gradient calculations');
	disp(toc);

	disp('Gradient mean of L1');
	disp(num2str(mean(L1)));

	disp('Gradient mean of L2');
	disp(num2str(mean(L2)));


end; % End of nargout if statement


%-----------------------------------------------------------------------------------------------------------------------------------------------%
% END OF FUNCTION
%-----------------------------------------------------------------------------------------------------------------------------------------------%
