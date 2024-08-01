%-----------------------------------------------------------------------------------------------%
% FUNCTION: Ax.m										%
% PURPOSE:  Calculate the left-hand side of the GIM equations used in the minimum		% 
% 	    residual approach.									%
% S. Miller, Oct. 9, 2018									%
%												%
%-----------------------------------------------------------------------------------------------%


function [ A ] = Ax(R, X, HX, Q, H, weights);

	% FUNCTION INPUTS:
	% R:		Model-data mismatch covariance matrix.
	% X:		Matrix of auxiliary variables.
	% HX:		Auxiliary variables passed through the forward model.
	% Q:		The Q covariance matrix (formatted as class kronMat)
	% H:		Class matvecH that points to the H matrices.
	% weights:	Current estimate for the kriging weights.	

	% FUNCTION OUTPUTS:
	% A:		Left-hand side of the GIM equations


%------------------------------------------------------%
% Separate out the weights into individual components  %
%------------------------------------------------------%

	% disp('-----------');
	disp('Compute A*x');
	tic;

        % Define the sizes of the different matrices
        p   = size(HX,2);
        n   = size(HX,1);	

        xi   = weights(1:(size(weights,1)-p));
        betas = weights((size(weights,1)-p+1):size(weights,1));


%----------------------------%
% Calculate the first line   %
%----------------------------%

        %! Note: Edit this section to match the actual format of H in your problem.

        % Calculate H^T * xi
	Hxi = H' * xi;

        % Calculate Q * H^T * xi
        % NOTE: This step is the slowest step for a year-long inverse model.
	% disp('Calculate Q * H^T * xi');
	QHxi = Q * Hxi;

	% Calculate H * Q * H' * xi
        % disp('Calculate HQHxi');
        % tic;
	HQHxi = H * QHxi;

	% Calculate R*xi
	% disp('Assemble the first line of Ax');
	Rxi = R*xi;

	% Calculate HXbeta
	HXbeta = HX*betas;
	
	% Sum the components together
	line1 = HQHxi + Rxi +  HXbeta;


%----------------------------%
% Calculate the second line  %
%----------------------------%

	% Calculate (HX)'*xi
	line2 = HX' * xi;


%--------------------------------------------%
% Finish calculations and return the result  %
%--------------------------------------------%

	A = [line1; line2];

	disp(toc);

%-----------------------------------------------------------------------------------------------%
% END OF FUNCTION										%
%-----------------------------------------------------------------------------------------------%

