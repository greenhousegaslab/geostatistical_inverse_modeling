%-------------------------------------------------------------------------------------------------------%
% FUNCTION: randeigdecomp										%
% PURPOSE: Compute the eigenvectors and eigenvalues of matrix A using a randomized approach.		%
% From Arvind Saibaba											%
%													%
%-------------------------------------------------------------------------------------------------------%


function [V,D] = randeigdecomp(Qchol, theta, n, Hpath, eigmax)
     
    % Implements a prototype randomized eigendecomposition
    %  
    % Inputs:
    % -------
    % 
    % Qchol   :     The Cholesky decomposition of Q (formatted as kronMat object).
    % theta   :     Matrix that contains covariance matrix parameters (see inversion_dual_eq.m) 
    % eigmax  :     target rank (i.e., number of eigenvalues/vectors to estimate)
    % H       :     H matrix (formatted as a matvecH class)
    % n       :     Total number of observations in the inverse problem.    

    % Outputs
    % -------- 
    % V       :     m      x eigmax       eigenvectors 
    % D       :     eigmax x eigmax       eigenvalues on the diagonals


	%------------------------------------------%
	% Determine the dimensions of the inputs   %
	%------------------------------------------%

    	m = size(Q,1);

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

	Y(:,l) = eigfun(Qchol, theta, n, H, Omega(:,l));

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
	temp(:,l) = eigfun(Qchol, theta, n, H, Q(:,l));
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

