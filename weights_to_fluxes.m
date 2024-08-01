%-----------------------------------------------------------------------------------------------%
% FUNCTION: weights_to_fluxes.m									%
% PURPOSE:  Calculate the fluxes given weights computed by the dual kriging equations.		%
% S. Miller, Nov. 22, 2018									%
%												%
%-----------------------------------------------------------------------------------------------%

%------------%
% NOTES:     %
%------------%


%---------------------%
% Begin function      %
%---------------------%

function [ shat ] = weights_to_fluxes(weights1,H,Q,X);

	% FUNCTION INPUTS:
	% weights1:	Estimated weights from the inverse model.
	% H: 		Object of class matvecH that points to the H matrix path.
	% Q: 		Q covariance matrix (of class kronMat).
	% X: 		Matrix of auxiliary variables.

	% FUNCTION OUTPUTS:
	% shat:		Posterior best estimate of the fluxes.


%------------------------------------%
% Estimate fluxes from the weights   %
%------------------------------------%

	disp('Estimate the fluxes from the weights');

        % Define the sizes of the different matrices
	p      = size(X,2);

	% Separate the weights into xi and beta
        xi   = weights1(1:(size(weights1,1)-p));
        betas = weights1((size(weights1,1)-p+1):size(weights1,1));;

        % Calculate H^T * xi
	Hxi = H'*xi;

        % Calculate Q * H^T * xi
	QHxi = Q * Hxi;

        % Create the estimate of the fluxes
        shat = X*betas + QHxi;

%-----------------------------------------------------------------------------------------------%
% END OF FUNCTION
%-----------------------------------------------------------------------------------------------%

