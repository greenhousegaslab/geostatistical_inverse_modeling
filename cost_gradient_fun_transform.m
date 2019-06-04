%-----------------------------------------------------------------------------------------------------------------------------------------------%
% FUNCTION; cost_gradient_fun_transform														%
% PURPOSE:  Compute the geostatistical inversion cost function. Use the Kronecker product to make these calculations more effeicient.		%
%	    This script will calculate both the cost function and the gradient (when the latter is requested).					%
% 	    This script uses a variable transformation to speed up convergence.									%
% S. Miller, Sept. 21, 2018															%
%																		%
%-----------------------------------------------------------------------------------------------------------------------------------------------%


	%-------------%
	% NOTES:      %
	%-------------%
		
	
function [ f,g ] = cost_gradient_fun_transform(Z, R, X, A, B, Dinv, Einv, CD, CE, Hpath, shat);

	% FUNCTION INPUTS:
	% Z		Observation vector
	% H		Footprints
	% R		Model-data mismatch matrix
	% X		Design matrix
	% Dinv		Inverse of the temporal correlation matrix (dimension nmonths x nmonths)
	% Einv		Inverse of the spatial correlation matrix (dimension land_grids x land_grids)
	% CD		Symmetric square root of D
	% CE		Symmetric square root of E
	% shat		Current guess for the fluxes 
	
	% FUNCTION OUTPUTS:
	% f		Value of the cost function
	% g		Value of the gradient
	
	
%--------------------------------------------%
% Set the dimensions of the inverse problem  %
%--------------------------------------------%

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

	% Calculate sqrt(Q)*shat
        Qx = [];

        for j = 1:ntimes;
        Qx1   = zeros(m1,1);
                for i = 1:ntimes;
                sel = (m1.*(i-1)+1):(i.*m1);
                Qx1 = Qx1 + shat(sel) .* CD(j,i);
                end; % End of i loop
        temp =  CE * Qx1;
        Qx   =  [Qx; temp];
        end; % End of i loop
        clear Qx1 temp;

	
	% Pass the flux estimate through the forward model
        %! Note: Edit this section to match the actual format of H in your problem.
	% disp('Calculate Hs');
        Hs = zeros(n,1);
        for j = 1:size(Dinv,1);
        load(strcat(Hpath,'H_',num2str(j),'.mat'));
        sel = (m1.*(j-1)+1):(j.*m1);
        Hs = Hs + H*Qx(sel,:);
        clear H;
        end;

	%  disp(toc);


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

	% Calculate  L1 = H' * (R \ (Z - Hs))] 
        % disp('Calculate H * (R \ (Z - Hs))');
        %! Note: Edit this section to match the actual format of H in your problem.
        L1 = [];
        for j = 1:size(Dinv,1);
        load(strcat(Hpath,'H_',num2str(j),'.mat'));
        L1 = [L1; H'* temp];
        end;

        % Multiply 'L1' by the symmetric square root of Q
	Qx = [];

        for j = 1:ntimes;
        Qx1   = zeros(m1,1);
                for i = 1:ntimes;
                sel = (m1.*(i-1)+1):(i.*m1);
                Qx1 = Qx1 + L1(sel,:) .* CD(j,i);
                end; % End of i loop
        temp =  CE * Qx1;
        Qx   =  [Qx; temp];
        end; % End of i loop
        clear Qx1 temp;
	L1 = Qx;


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
