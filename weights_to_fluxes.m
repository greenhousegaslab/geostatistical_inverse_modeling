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

function [ shat ] = weights_to_fluxes(weights1,Hpath,D,E,X);

	% FUNCTION INPUTS:
	% weights1:	Estimated weights from the inverse model.
	% Hpath: 	Path to the H matrix.
	% D: 		Matrix of temporal covariances.
	% E:		Matrix of spatial covariances.
	% X: 		Matrix of auxiliary variables.

	% FUNCTION OUTPUTS:
	% shat:		Posterior best estimate of the fluxes.


%------------------------------------%
% Estimate fluxes from the weights   %
%------------------------------------%

	disp('Estimate the fluxes from the weights');

        % Define the sizes of the different matrices
        m1     = size(E,1);
        ntimes = size(D,1);
	m      = m1.*ntimes;
	p      = size(X,2);

	% Separate the weights into xi and beta
        xi   = weights1(1:(size(weights1,1)-p));
        betas = weights1((size(weights1,1)-p+1):size(weights1,1));;

        % Calculate H * xi
        Hxi = [];
        for j = 1:size(D,1);
        H = load(strcat(Hpath,'H_',num2str(j),'.mat'));
       	H = H.H; 
	Hxi = [Hxi; H'*xi];
        end;

        % Calculate Q * H^T * xi
        QHxi = [];

        for j = 1:ntimes;
        A1 = zeros(m1,1);
                for i = 1:size(D,1);
                sel = (m1.*(i-1)+1):(i.*m1);
                A1 = A1 + Hxi(sel,1) .* D(j,i);
                end; % End of i loop
        temp = E * A1;
        QHxi = [QHxi; temp];
        end; % End of j loop
        A1 = []; temp = [];

        % Create the estimate of the fluxes
        shat = X*betas + QHxi;


%-----------------------------------------------------------------------------------------------%
% END OF FUNCTION
%-----------------------------------------------------------------------------------------------%

