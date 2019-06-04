%-------------------------------------------------------------------------------------------------------%
% FUNCTION: randeigdecomp										%
% PURPOSE: Compute the eigenvectors and eigenvalues of matrix A using a randomized approach.		%
% From Arvind Saibaba											%
%													%
%-------------------------------------------------------------------------------------------------------%


function [V,D] = randeigdecomp(CD, CE, theta, n, Hpath, eigmax)
     
    % Implements a prototype randomized eigendecomposition
    %  
    % Inputs:
    % -------
    % 
    % CD      :     The Cholesky decomposition of matrix D.
    % CE      :     The Cholesky decomposition of matrix E.
    % theta   :     Matrix that contains covariance matrix parameters (see inversion_dual_eq.m) 
    % eigmax  :     target rank (i.e., number of eigenvalues/vectors to estimate)
    % Hpath   :     Path to the H matrix (e.g., sensitivity matrix or matrix of footprints)
    % n       :     Total number of observations in the inverse problem.    

    % Outputs
    % -------- 
    % V       :     m      x eigmax       eigenvectors 
    % D       :     eigmax x eigmax       eigenvalues on the diagonals


	%------------------------------------------%
	% Determine the dimensions of the inputs   %
	%------------------------------------------%

    	m = size(CD,1).*size(CE,1);

	% Set the oversampling parameter
    	p = 20;
    	% p = 200;


	%-----------------------------------------------------------%
	% Create a random matrix and multiply it by input matrix A  %
	%-----------------------------------------------------------%

    	Omega = randn(m,eigmax+p);
    	% Y     = A*Omega;
    

	%------------------------%
	% Multiply A and Omega   %
	%------------------------%

	Y = zeros(m,eigmax+p);

	parfor l = 1:(eigmax+p);

	Y(:,l) = eigfun(CD, CE, theta, n, Hpath, Omega(:,l));

	end;


	%------------------------------------------------%
	% Create QR decomposition of Y and mutiply by A  %
	%------------------------------------------------%

	% Y has dimensions m x (eigmax+p)
	% Q is small: dimensions (eigmax+p) x (eigmax+p)
	% 0 option in qr function: "If m > n, only the first n columns of Q and the first n rows of R are computed"
    	[Q,~] = qr(Y, 0); 
	clear Y;

	% Option 1: Requires computations with A. Is exact
    	%% T     = Q'*(A*Q);
	temp = zeros(m,eigmax+p);
	parfor l = 1:(eigmax + p);
	temp(:,l) = eigfun(CD, CE, theta, n, Hpath, Q(:,l));
	end;
	T = Q' * temp;
	clear temp;

	% Option 2: Approximate. Does not require computations with A
	% This approach is really inaccurate for large inverse problems.
	% The estimated eigenvalues are always too small (in my experience)
	% T = (Q'*Y)*inv(Q'*Omega);
    

	%---------------------------------------------%
	% Find the eigenvectors of T and order them   %
	%---------------------------------------------%
 
   	[U,D]   = eig(T,'vector');
    	[D,ind] = sort(D, 'descend');
    

	%--------------------------------------------------%
	% Multiply Q by the eigenvectors calculated above  %
	%--------------------------------------------------%

	V = Q*U(:,ind);
	clear Q U;

	%------------------------------------------------%
	% Find the eigenvectors (V) and eigenvalues (D)  %
	%------------------------------------------------%
    
    	V = V(:,1:eigmax);
    	D = diag(D(1:eigmax));


%-------------------------------------------------------------------------------------------------------%
end
%-------------------------------------------------------------------------------------------------------%

