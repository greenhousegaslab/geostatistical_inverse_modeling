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
		
	
function [ f,g ] = cost_gradient_fun(Z, R, X, B, Dinv, Einv, Hpath, shat);

	% FUNCTION INPUTS:
	% Z		Observation vector
	% H		Footprints
	% R		Model-data mismatch matrix
	% X		Design matrix
	% Dinv		Inverse of the temporal correlation matrix (D)
	% Einv		Inverse of the spatial correlation matrix (E)

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
        Hs = zeros(n,1);
        for j = 1:size(Dinv,1);
        load(strcat(Hpath,'H_',num2str(j),'.mat'));
        sel = (m1.*(j-1)+1):(j.*m1);
        Hs = Hs + H*shat(sel,:);
        clear H;
        end;


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
	A = [];
	
	for j = 1:size(Dinv,1);
	A1 = zeros(m1,1);
		for i = 1:size(Dinv,1);
		sel = (m1.*(i-1)+1):(i.*m1);
		A1 = A1 + shat(sel) .* Dinv(j,i);	
		end; % End of i loop
	temp = Einv * A1;
	A = [A; temp];
	end; % End of i loop
	clear A1 temp;

	
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
        L1 = [];
        for j = 1:size(Dinv,1);
        load(strcat(Hpath,'H_',num2str(j),'.mat'));
        L1 = [L1; H'*temp];
        end;


%------------------------------------------------%
% Calculate the prior component of the gradient  %
%------------------------------------------------%

	% Calculate G * shat
	% G * shat = A - B * inv( X' * B) * B' * shat
	
	%---------------%
	% Calculate A   %
	%---------------%
	
	% A = inv(Q) * shat
	A = [];
	
	for j = 1:size(Dinv,1);
	A1 = zeros(m1,1);
		for i = 1:size(Dinv,1);
		sel = (m1.*(i-1)+1):(i.*m1);
		A1 = A1 + shat(sel) .* Dinv(j,i);	
		end; % End of i loop
	temp = Einv * A1;
	A = [A; temp];
	end; % End of i loop
	clear A1 temp;
	
	
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
