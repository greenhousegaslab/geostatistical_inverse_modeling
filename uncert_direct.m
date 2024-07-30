%---------------------------------------------------------------------------------------%
% FUNCTION: uncert_direct.m								%
% PURPOSE:  Calculate the uncertainties using direct calculations.			%
% S. Miller, Nov. 28, 2018								%
%	Updated July 24, 2024								%
%											%
%---------------------------------------------------------------------------------------%


%----------%
% NOTES:   %
%----------%

        % Note: In this script, we assume that the fluxes are being estimated at a 3-hourly time resolution
        % and have units of micromol m-2 s-1. Furthermore, we assume that the results are being reported
        % in units of petagrams carbon. If that is not the case, make sure to edit the units conversion
        % where noted in the script below.


%------------------%
% Begin function   %
%------------------%

function [ uncert ] = uncert_direct(psi,HX,X,Q,H,selu);

	% FUNCTION INPUTS:
	% psi:		HQH' + R
	% HX:		The sensitivity matrix (H) multiplied by the auxiliary variables (X)
	% X:		Matrix of auxiliary variables.
	% Q:		The Q covariance matrix. This object should be formatted using the kronMat class.
	% H: 		An object of class "matvecH" that points to the H matrix files.
	% selu: 	This vector contains the areas of each model grid box.
	%		The elements of this vector should be set to 0 for all model grid boxes that 
	%		should be excluded from the uncertainty calculations.

	% OUTPUTS:
	% uncert: 	95% confidence interval = posterior best estimate +/- this number 

        % Equation for the uncertainties: 
        % See the appendix of Kitanidis 1996 for the full equation
        % https://link.springer.com/content/pdf/10.1007%2FBF01581870.pdf
        % V_post = Q - (HQ)'*P*(HQ) - X*B*X' - A'*HQ - (HQ)'*A
        % Naming: V_post = V1 - V2 - V3 - V4 - V5 
        % inv([psi HX; HX' 0]) = [P A; A' B]


%------------------------------------------------%
% Set several parameters required by the script  %
%------------------------------------------------%

        p   = size(X,2);
        n   = size(HX,1);
        m1  = size(selu,1);
        ntimes = size(Q,1)./m1;


%----------------------------------%
% Create unit conversion vectors   %
%----------------------------------%

        % Repeat the unit conversion vector based upon the number of time periods in the inverse model
        agvec = repmat(selu,ntimes,1);


%--------------------------%
% Calculate P, A, and B    %
%--------------------------%

        % Create P, A, and B
        temp = [psi HX; HX' zeros(p,p)];
        clear psi HX;
        temp = temp \ speye(n+p,n+p);
        P = temp(1:n,1:n);
        A = temp(1:n,(n+1):(n+p));
        B = temp((n+1):(n+p),(n+1):(n+p));
        clear temp;
	% save(strcat(outpath,'invpsi.mat'),'P','A','B','-v7.3');


%---------------%
% Calculate V1  %
%---------------%

        % Multiply Q by aggregation vector
	Qag = Q * agvec; 

        % Finish calculating (agvec)' * Q * agvec
        V1 = agvec' * Qag;


%----------------%
% Calculate V2   %
%----------------%

        % H*Qag
	HQx = H * Qag;

        % Multiply by P
        PHQx = P*HQx;
        clear HQx;

        % Multiply by H'
	HPHQx = H' * PHQx;

        % Mulitiply by Q' (Same as multiplying by Q since Q and Q' are symmetric)
	QHPHQx = Q*HPHQx;

        % Multiply by agvec'
        V2 = agvec' * QHPHQx;


%----------------%
% Calculate V3   %
%----------------%

        % Multiply X'*agvec
        Xa = X'*agvec;

        % Multiply by B
        BXa = B*Xa;
        clear Xa;

        % Multiply by X
        XBXa = X*BXa;
        clear BXa;

        % Multiply by the aggregation vector
        V3 = agvec' * XBXa;


%----------------%
% Calculate V4   %
%----------------%

        % X*A'*HQ   (mxp)*(pxn)*(nxm)

	% Multiply H by (Q*agvec)
	HQx = H * Qag; 

        % Multiply by X*A'
        XAQag = X*(A'*HQx);

        % Multiply by agvec
        V4 = agvec'*XAQag;


%----------------%
% Calculate V5   %
%----------------%

        % Q * H'*A*X' (mxn)*(nxp)*(pxm)

        % A*X'*agvec
        AXa = A * (X' * agvec);

        % Multiply by H'
	HAXa = H' * AXa;

        % Multiply by Q
	QHAXa = Q * HAXa;

        % Multiply by agvec' 
        V5 = agvec' * QHAXa;


%----------------------------------%
% Put all of the pieces together   %
%----------------------------------%

        uncert = V1 - V2 - V3 - V4 - V5;

        uncert = sqrt(uncert);

        % Do final unit conversions
        % Convert from micromol to Pg C 
        % Convert micromol to mol
        uncert = uncert ./ 1e6;

        % Convert mol to gram C
        uncert = uncert .* 12;

        % Convert gram to petagram
        uncert = uncert .* 1e-15;

        % Convert to 95% confidence interval
        uncert = 1.96 .* uncert;

        disp('95% confidence interval =');
	disp('posterior best estimate +/- this number:');
        disp(num2str(uncert));


%---------------------------------------------------------------------------------------%
% END OF FUNCTION									%
%---------------------------------------------------------------------------------------%

