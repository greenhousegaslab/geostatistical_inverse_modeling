%-----------------------------------------------------------------------------------------------------------------------%
% FUNCTION: HQHR.m													%
% PURPOSE: Calculate (HQH' + R) and return the result. 									%
% S. Miller, July 13, 2015												%
%															%
%-----------------------------------------------------------------------------------------------------------------------%


	%-----------%
	% Notes:    %
	%-----------%
		

function [ psi ] = HQHR(R, D, E, Hpath);

	% INPUTS:
	% H:		Sensitivity matrix (dimensions n x m)
	% R:		Model-data mismatch covariance matrix.
	% D: 		Matrix that describes temporal covariances in the fluxes
	%		(This matrix should be multiplied by the variance in Q.)
	% E: 		Matrix that describes spatial covariances in the fluxes
	% Hpath:	Path to the H matrix. See the readme file for more info on 
	%		options for handling H.

	% OUTPUTS:
	% psi:		(HQH' + R)


%-----------------------------%
% Set initial size parameters %
%-----------------------------%

	n       = size(R,1);
	m       = size(D,1) .* size(E,1);
	m1      = size(E,1);
        ntimes  = size(D,1);


%---------------------%
% Calculate HQH       %
%---------------------%

	disp('Calculate HQH');

	psi = sparse(n,n);

	for j = 1:ntimes;
	disp(j);
	for i = 1:ntimes;
	if D(i,j) > 0 | i==1;
        	load(strcat(Hpath,'H_',num2str(i),'.mat'));
		if i==1; 
		HQ1 = H .* D(i,j);
		else;
		HQ1 = HQ1 + H .* D(i,j);
		end; % End of if statement
	end; % End of if statement
	end; % End of i loop

	HQ1 = HQ1 * E;

        load(strcat(Hpath,'H_',num2str(j),'.mat'));
	HQ1 = HQ1 * H';

	psi = psi + HQ1;

	end; % End of j loop
	clear HQ1;
	

%---------------------------%
% Add R to complete psi     %
%---------------------------%

	disp('Add R to HQH');
	psi = psi + R;


%-----------------------------------------------------------------------------------------------------------------------%
% END OF FUNCTION													%
%-----------------------------------------------------------------------------------------------------------------------%
