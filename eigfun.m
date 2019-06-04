%-----------------------------------------------------------------------------------------------%
% FUNCTION: eigfun.m										%
% PURPOSE: Function for computing Ax where A is the matrix we want eigenvectors/values for.	%
% S. Miller, Nov. 15, 2018									%
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

function [ Ax ] = eigfun(CD, CE, theta, n, Hpath, x1);

	% FUNCTION INPUTS:
	% theta: 	First value of this vector should be sigma_R
	% n: 		Number of observations
	% CD: 		Cholesky decomposition of the D matrix
	% CE: 		Cholesky decomposition of the E matrix
	% Hpath: 	Path to the footprint matrices
	% x1: 		Unknown vector that is provided by the "eigs" function in Matlab.
	%     		or by the random eigenvector decomposition function (randeigdecomp.m)

disp('Computing matrix-vector products for eigenvector/value calculation');
tic;

%--------------------------------%
% Define the problem dimensions  %
%--------------------------------%

	ntimes = size(CD,1);
	m = size(CD,1).*size(CE,1);
	m1 = m ./ ntimes; % Divide m by the total number of time periods in the inversion


%--------------------%
% Compute chol(Q)*x  %
%--------------------%

	% I've included a transpose here.
	CDt = CD';
	Qx = [];

	for j = 1:ntimes;
	Qx1   = zeros(m1,1);
		for i = 1:ntimes;
		sel = (m1.*(i-1)+1):(i.*m1);
		Qx1 = Qx1 + x1(sel,:) .* CDt(j,i);
		end; % End of i loop
	temp = 	CE' * Qx1;
	Qx   =  [Qx; temp];
	end; % End of i loop
	clear Qx1 temp;


%-----------------%
% Multiply by H   %
%-----------------%

        %! Note: Edit this section to match the actual format of H in your problem.

        HQx = zeros(n,1);
        for j = 1:ntimes;
        load(strcat(Hpath,'H_',num2str(j),'.mat'));
        sel = (m1.*(j-1)+1):(j.*m1);
        HQx = HQx + H*Qx(sel,:);
        clear H;
        end;
	clear Qx1;


%--------------------------%
% Multiply by R^(-1)       %
%--------------------------%

	RHQx = (1./(theta(1).^2)) .* HQx;
	clear HQx;
	

%-------------------%
% Multiply by H^T.  %
%-------------------%

        %! Note: Edit this section to match the actual format of H in your problem.

	HRHQx = [];
	for j = 1:ntimes;
	load(strcat(Hpath,'H_',num2str(j),'.mat'));
	HRHQx = [ HRHQx; H'* RHQx ];
	end;
	clear RHQx;


%-----------------------%
% Multiply by chol(Q)   %
%-----------------------%

	Ax = [];

	for j = 1:ntimes;
	Ax1   = zeros(m1,1);
		for i = 1:ntimes;
		sel = (m1.*(i-1)+1):(i.*m1);
		Ax1 = Ax1 + HRHQx(sel) .* CD(j,i);
		end; % End of i loop
	temp = 	CE * Ax1;
	Ax   =  [Ax; temp];
	end; % End of i loop
	clear Ax1 temp;

	disp(toc);


%-----------------------------------------------------------------------------------------------%
% END OF FUNCTION
%-----------------------------------------------------------------------------------------------%

