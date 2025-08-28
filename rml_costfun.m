%---------------------------------------------------------------------------------------------------------------------------------------%
% FUNCTION: rml_costfun.m														%
% PURPOSE: Compute the value of the RML cost function											%
% S. Miller, July 27, 2015														%
%	Updated Aug. 28, 2025														%
%																	%
%---------------------------------------------------------------------------------------------------------------------------------------%

%-----------%
% Notes:    %
%-----------%

	% See rml.m for more details (and to run RML).


function [ L ] = rml_costfun(theta, Z, X, H, deltamat, days);

        disp('Evaluate RML cost function');

        % Keep track of the amount of time it takes to evaluate the cost function.
        tic;


%------------------%
% Initial setup    %
%------------------%

        disp('Covariance parameters:');
        disp(theta);

	n 	= length(Z);
	m 	= size(X,1);
        p       = size(X,2);


%--------------------------------%
% Create the D and E matrices    %
%--------------------------------%

        disp('Create the E and D matrices');

        % The E and D matrices are used to construct the Q covariance matrix. 
        % The E matrix describes how covariances in Q decay in space, and the D matrix describes how covariances in Q decay in time.
        % See the following paper by Vineet Yadav for more details on the E and D matrices: http://www.geosci-model-dev.net/6/583/2013/gmd-6-583-2013.html
        % The E matrix has dimensions (r x r)
        % The D matrix has dimensions (q x q)

        %------------%
        % Create E   %
        %------------%

        % No need to edit this section
        % Note: I use a spherical covariance model because of its sparse properties.
	% Comment out the exponential model below if you would rather use that model.

        % Multiply the distance matrix by the covariance model
        % Exponential covariance model
        % decaymat = exp(-1 .* deltamat ./ theta(3));

        % Spherical covariance model
        E                      = 1 - 1.5 .* (deltamat ./theta(3))  + 0.5 .* (deltamat.^3 ./ theta(3).^3);
        E(deltamat > theta(3)) = 0;


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


%------------------------%
% Create the R matrix    %
%------------------------%

        disp('Create R');

        % No need to edit this section, unless you want a more complicated setup for R.
        R = (theta(1).^2) .* ones(n,1);
        R = spdiags(R,0,n,n);


%--------------------%
% Compute psi        %
%--------------------%

	D = (sigmaQ*sigmaQ') .* D;
        psi = HQHR(H, R, D, E);


%---------------------------------------%
% Compute term 1 of the cost function   %
%---------------------------------------%

        L1 = 0.5.*2.*sum(log(diag(chol(psi))));


%--------------------------------%
% Calculate psi\H and psi\z      %
%--------------------------------%

	% NOTE: These two calculations are really the time bottle-neck in this script!
	%	I don't think there's any way to avoid these calculations, though.
	HX = H*X;
        psiHX = psi\HX;
        psiZ = psi\Z;
        clear psi;


%---------------------------------------%
% Compute term 2 of the cost function   %
%---------------------------------------%

        % Here, we compute the determinant of a relatively small matrix. 
        % Hence, we don't need to use any clever means to calculate the determinant

        XHpsiHX = (HX)' * psiHX;
	L2 = 0.5.*log(det(XHpsiHX));
	clear HX;


%---------------------------------------%
% Compute term 3 of the cost function   %
%---------------------------------------%

        L3 = 0.5 .* (Z'*psiZ);
        ZpsiHX = (Z' * psiHX);
	L4 = -0.5 .* ZpsiHX * (XHpsiHX \ (ZpsiHX'));


%-------------------------------------%
% Finish computing the cost function  %
%-------------------------------------%

        L = L1 + L2 + L3 + L4;

        disp('Cost function evaluated. Value:');
        disp(num2str(round(L)));

	disp('Time to compute the cost function (sec)');
	disp(num2str(round(toc)));


%---------------------------------------------------------------------------------------------------------------------------------------%
% END OF FUNCTION															%
%---------------------------------------------------------------------------------------------------------------------------------------%
