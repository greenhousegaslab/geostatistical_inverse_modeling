%-----------------------------------------------------------------------------------------------------------------------%
% FUNCTION: HQHR.m													%
% PURPOSE: Calculate (HQH' + R) and return the result. 									%
% S. Miller, July 13, 2015												%
%															%
%-----------------------------------------------------------------------------------------------------------------------%


	%-----------%
	% Notes:    %
	%-----------%
		

function [ psi ] = HQHR(R, D, E, H);

	% INPUTS:
	% R:		Model-data mismatch covariance matrix.
	% D: 		Matrix that describes temporal covariances in the fluxes
	%		(This matrix should be multiplied by the variance in Q.)
	% E: 		Matrix that describes spatial covariances in the fluxes
	% H:		The H matrix (i.e., footprints). We recommend using the function "matvecH" to create this object.

	% OUTPUTS:
	% psi:		(HQH' + R)


%-----------------------------%
% Set initial size parameters %
%-----------------------------%

	n       = size(R,1);
	m       = size(D,1) .* size(E,1);
	m1      = size(E,1);
        ntimes  = size(D,1);

%--------------------------------------%
% Create Q using the kronMat function. %
%--------------------------------------%

	Q = kronMat(D,E);


%---------------------%
% Calculate HQH       %
%---------------------%

	disp('Calculate HQH');

	psi = (H*Q)';
	psi = H * psi;
	

%---------------------------%
% Add R to complete psi     %
%---------------------------%

	disp('Add R to HQH');
	psi = psi + R;


%-----------------------------------------------------------------------------------------------------------------------%
% END OF FUNCTION													%
%-----------------------------------------------------------------------------------------------------------------------%
