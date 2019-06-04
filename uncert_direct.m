%---------------------------------------------------------------------------------------%
% FUNCTION: uncert_direct.m								%
% PURPOSE:  Calculate the uncertainties using direct calculations.			%
% S. Miller, Nov. 28, 2018								%
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

function [ uncert ] = uncert_analytical(psi,HX,D,E,X,selu,Hpath);

	% FUNCTION INPUTS:
	% psi:		HQH' + R
	% HX:		The sensitivity matrix (H) multiplied by the auxiliary variables (X)
	% D:		Matrix that describes temporal covariance in the fluxes.
	%		(Note that D should be multiplied by the variance in Q.)
	% E:		Matrix that describes spatial covariance in the fluxes.
	% X:		Matrix of auxiliary variables.
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
        m1  = size(E,1);
        m   = size(E,1).*size(D,1);
	ntimes = size(D,1);


%----------------------------------%
% Create unit conversion vectors   %
%----------------------------------%

        % Convert time from (1/s) to (1/time period) [ in this case, 3 hours ]
        selu = selu.*60.*60.*3;

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
        Qag = [];

        for j = 1:ntimes;
        A1 = zeros(m1,1);
                for i = 1:ntimes;
                sel = (m1.*(i-1)+1):(i.*m1);
                A1 = A1 + agvec(sel) .* D(j,i);
                end; % End of i loop
        temp = E * A1;
        Qag = [Qag; temp];
        end; % End of i loop
        clear A1 temp;

        % Finish calculating (agvec)' * Q * agvec
        V1 = agvec' * Qag;


%----------------%
% Calculate V2   %
%----------------%

        %! Note: Edit this section to match the actual format of H in your problem.

        % H*Qag
        HQx = zeros(n,1);
        for j = 1:ntimes;
        load(strcat(Hpath,'H_',num2str(j),'.mat'));
        sel = (m1.*(j-1)+1):(j.*m1);
        HQx = HQx + H*Qag(sel,:);
        clear H;
        end;
        clear Qx1;

        % Multiply by P
        PHQx = P*HQx;
        clear HQx;

        % Multiply by H'
        HPHQx = [];
        for j = 1:ntimes;
	load(strcat(Hpath,'H_',num2str(j),'.mat'));
	HPHQx = [ HPHQx; H'* PHQx ];
        end;
        clear PHQx;

        % Mulitiply by Q' (Same as multiplying by Q)
        Qag = [];

        for j = 1:ntimes;
        A1 = zeros(m1,1);
                for i = 1:ntimes;
                sel = (m1.*(i-1)+1):(i.*m1);
                A1 = A1 + HPHQx(sel,:) .* D(j,i);
                end; % End of i loop
        temp = E * A1;
        Qag = [Qag; temp];
        end; % End of i loop
        clear A1 temp;

        % Multiply by agvec'
        V2 = agvec' * Qag;


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

        % Multiply Q*agvec
        Qag = [];

        for j = 1:ntimes;
        A1 = zeros(m1,1);
                for i = 1:ntimes;
                sel = (m1.*(i-1)+1):(i.*m1);
                A1 = A1 + agvec(sel) .* D(j,i);
                end; % End of i loop
        temp = E * A1;
        Qag = [Qag; temp];
        end; % End of i loop
        clear A1 temp;

	% Multiply by H
        HQx = zeros(n,1);
        for j = 1:ntimes;
        load(strcat(Hpath,'H_',num2str(j),'.mat'));
	HQx = HQx + H*Qag(sel,:);
        clear H;
        end;
        clear Qx1;

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
        HAXa = [];
        for j = 1:ntimes;
        load(strcat(Hpath,'H_',num2str(j),'.mat'));
	HAXa = [ HAXa; H'* AXa ];
        end;
        clear AXa;

        % Multiply by Q
        Qag = [];

        for j = 1:ntimes;
        A1 = zeros(m1,1);
                for i = 1:ntimes;
                sel = (m1.*(i-1)+1):(i.*m1);
                A1 = A1 + HAXa(sel,:) .* D(j,i);
                end; % End of i loop
        temp = E * A1;
        Qag = [Qag; temp];
        end; % End of i loop
        clear A1 temp;

        % Multiply by agvec' 
        V5 = agvec' * Qag;


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

