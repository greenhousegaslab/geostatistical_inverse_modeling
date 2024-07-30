%-----------------------------------------------------------------------------------------------------------------------------------------------%
% FUNCTION; cost_gradient_fun_transform														%
% PURPOSE:  Compute the geostatistical inversion cost function. Use the Kronecker product to make these calculations more effeicient.		%
%	    This script will calculate both the cost function and the gradient (when the latter is requested).					%
% 	    This script uses a variable transformation to speed up convergence.									%
% S. Miller, Sept. 21, 2018															%
%	Updated July 24, 2024															%
%																		%
%-----------------------------------------------------------------------------------------------------------------------------------------------%


	%-------------%
	% NOTES:      %
	%-------------%
		
	
function [ f,g ] = cost_gradient_fun_transform(Z, R, X, A, B, Qinv, Qsqrt, H, shat);

	% FUNCTION INPUTS:
	% Z		Observation vector
	% H		Footprints
	% R		Model-data mismatch matrix
	% X		Design matrix
	% Qinv		Inverse of the Q covariance matrix (formatted as a kronMat class)
	% Qsqrt		Matrix square root of Q (formatted as a kronMat class)
	% shat		Current guess for the fluxes 
	
	% FUNCTION OUTPUTS:
	% f		Value of the cost function
	% g		Value of the gradient
	

%---------------------------------------%
%***** CALCULATE THE COST FUNCTION *****%
%---------------------------------------%

disp('Calculating the cost function');


%------------------------------------------------------------%
% Pass the current flux estimate through the forward model   %
%------------------------------------------------------------%

	tic;

	% Calculate sqrt(Q)*shat
	Qx = Qsqrt * shat;
	
	% Pass the flux estimate through the forward model
	% disp('Calculate Hs');
	Hs = H * Qx;


%-------------------------------------------%
% Calculate the likelihood/data component   %
%-------------------------------------------%

	% Fastest formulation
	L1 = sum((Z - Hs).^2./diag(R));


%---------------------------------%
% Calculate the prior component   %
%---------------------------------%
	
	% Detailed calculations:
	% sQinv = inv(sqrtm(Q));
	% L2 = shat'* [I + (sQinv * X) * ((X' * Qinv * X) \ (sQinv * X)') ] * shat
	
	Gshat = shat - A * ((X' * B) \ (A' * shat));	

	L2 = shat' * Gshat;
	clear Gshat;


%-------------------------------------------------%
% Add up the two components of the cost function  %
%-------------------------------------------------%

	f = 0.5 .* (L1 + L2);

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

	tic;
	temp = (R \ (Z - Hs));

	% Calculate  L1 = H' * (R \ (Z - Hs))] 
        % disp('Calculate H * (R \ (Z - Hs))');
        %! Note: Edit this section to match the actual format of H in your problem.
	L1 = H' * temp;

        % Multiply 'L1' by the symmetric square root of Q
	L1 = Qsqrt * L1;


%------------------------------------------------%
% Calculate the prior component of the gradient  %
%------------------------------------------------%

	% Calculate G * shat
        % L2 = [I + (sQinv * X) * ((X' * Qinv * X) \ (sQinv * X)') ] * shat

       L2 = shat - A * ((X' * B) \ (A' * shat));	


%--------------------------------------------%
% Add up the two components of the gradient  %
%--------------------------------------------%

	g = L2 - L1;

	disp('Mean of the gradient');
	disp(num2str(mean(g)));

	disp('Time elapsed for gradient calculations');
	disp(toc);


end; % End of nargout if statement


%-----------------------------------------------------------------------------------------------------------------------------------------------%
% END OF FUNCTION
%-----------------------------------------------------------------------------------------------------------------------------------------------%
