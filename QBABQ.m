%-----------------------------------------------------------------------------------------------%
% FUNCTION: QBABQ.m										%
% PURPOSE: Calculate Q^(1/2)*B*A*B'*Q^(1/2)*x2 where x2 is some vector or matrix. 		%
% S. Miller, Nov. 20, 2018									%
%												%
%-----------------------------------------------------------------------------------------------%


function [ QBABQx ] = QBABQ(Qchol, eigvec, eigval, x2);

	% FUNCTION INPUTS:
	% Qchol:	Cholesky decomposition of Q (formatted as a kronMat object)
	% eigvec:	Eigenvectors of the prior preconditioned Hessian matrix	
	% eigval:	Eigenvalues of the prior preconditioned Hessian matrix
	% x2:		Vector that gets multiplied by Q^(1/2)*B*A*B'*Q^(1/2)

	% FUNCTION OUTPUTS:
	% QBABQx: Q^(1/2)*B*A*B'*Q^(1/2)*x2


% Calculate Q^(1/2)*x2
        Qx2 = Qchol * x2;

% Multiply by B*A*B'
        BABQx = (eigvec*(eigval*(eigvec'*Qx2)));
        clear Qx2;

% Multiply by Q^(1/2)
        QBABQx = Qchol * BABQx;


%-----------------------------------------------------------------------------------------------%
% END OF FUNCTION
%-----------------------------------------------------------------------------------------------%


