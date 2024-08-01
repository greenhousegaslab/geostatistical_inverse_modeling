%-----------------------------------------------------------------------------------------------%
% FUNCTION: eigfun.m										%
% PURPOSE: Function for computing Ax where A is the matrix we want eigenvectors/values for.	%
% S. Miller, Nov. 15, 2018									%
%	Updated July 24, 2024									%
%												%
%-----------------------------------------------------------------------------------------------%

%-------------------%
% Function notes:   %
%-------------------%

	% This function will calculate chol(Q)*H^T*R^(-1)*H*chol(Q)*x
	% where x is some vector specified by the "eigs" function in Matlab

	% This function requires passing a vector through the forward atmospheric model 
	% and a vector through the adjoint model (H'). In this script, I've assumed that H
	% can be read in and stored in computer memory. For many large problems, that is not the case.
	% Rather, it may be advisable to read in the H matrix in pieces or strips (if there is an explicit H matrix).
	% Alternately, one can pass a vector through the forward model and then pass a vector through the adjoint
	% model without explicitely formulating H. These sections where these calculations are required are
	% denoted below.


%------------------%
% Begin function   %
%------------------%

function [ Ax ] = eigfun(Qchol, theta, n, H, x1);

	% FUNCTION INPUTS:
	% Qchol:	Cholesky decomposition of Q (formatted as a kronMat object)
	% theta: 	First value of this vector should be sigma_R
	% n: 		Number of observations
	% H: 		H as a matvecH object
	% x1: 		Unknown vector that is provided by the "eigs" function in Matlab.
	%     		or by the random eigenvector decomposition function (randeigdecomp.m)

disp('Computing matrix-vector products for eigenvector/value calculation');

%--------------------%
% Compute chol(Q)*x  %
%--------------------%

	Qx = Qchol * x1;


%-----------------%
% Multiply by H   %
%-----------------%

        HQx = H * Qx;


%--------------------------%
% Multiply by R^(-1)       %
%--------------------------%

	RHQx = (1./(theta(1).^2)) .* HQx;
	clear HQx;
	

%-------------------%
% Multiply by H^T.  %
%-------------------%

        %! Note: Edit this section to match the actual format of H in your problem.

	HRHQx = H' * RHQx;


%-----------------------%
% Multiply by chol(Q)   %
%-----------------------%

	Ax = Qchol * HRHQx;


%-----------------------------------------------------------------------------------------------%
% END OF FUNCTION
%-----------------------------------------------------------------------------------------------%

