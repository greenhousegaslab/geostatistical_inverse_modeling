%-----------------------------------------------------------------------------------------------%
% FUNCTION: QBABQ.m										%
% PURPOSE: Calculate Q^(1/2)*B*A*B'*Q^(1/2)*x2 where x2 is some vector or matrix. 		%
% S. Miller, Nov. 20, 2018									%
%												%
%-----------------------------------------------------------------------------------------------%


function [ QBABQx ] = QBABQ(CD, CE, eigvec, eigval, x2);

	% FUNCTION INPUTS:
	% CD:		Cholesky decomposition of D (matrix that describes temporal covariances)
	% CE:		Cholesky decomposition of E (matrix that describes spatial covariances)
	% eigvec:	Eigenvectors of the prior preconditioned Hessian matrix	
	% eigval:	Eigenvalues of the prior preconditioned Hessian matrix
	% x2:		Vector that gets multiplied by Q^(1/2)*B*A*B'*Q^(1/2)

	% FUNCTION OUTPUTS:
	% QBABQx: Q^(1/2)*B*A*B'*Q^(1/2)*x2


% Important parameters for the function
        ntimes = size(CD,1);
        m = size(CD,1).*size(CE,1);
        m1 = m ./ ntimes; % Divide m by the total number of time periods in the inversion

% Calculate Q^(1/2)*x2
        Qx2 = [];

        for j = 1:ntimes;
        Qx1   = zeros(m1,1);
                for i = 1:ntimes;
                sel = (m1.*(i-1)+1):(i.*m1);
                Qx1 = Qx1 + x2(sel,:) .* CD(j,i);
                end; % End of i loop
        temp =  CE * Qx1;
        Qx2   =  [Qx2; temp];
        end; % End of i loop
        clear Qx1 temp;

% Multiply by B*A*B'
        BABQx = (eigvec*(eigval*(eigvec'*Qx2)));
        clear Qx2;

% Multiply by Q^(1/2)
        QBABQx = [];
	CDt = CD';

        for j = 1:ntimes;
        Qx1   = zeros(m1,1);
                for i = 1:ntimes;
                sel = (m1.*(i-1)+1):(i.*m1);
                Qx1 = Qx1 + BABQx(sel,:) .* CDt(j,i);
                end; % End of i loop
        temp =  CE' * Qx1;
        QBABQx   =  [QBABQx; temp];
        end; % End of i loop
        clear Qx1 temp BABQx;


%-----------------------------------------------------------------------------------------------%
% END OF FUNCTION
%-----------------------------------------------------------------------------------------------%


