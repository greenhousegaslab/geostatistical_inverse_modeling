%-----------------------------------------------------------------------------------------------%
% FUNCTION: Ax.m										%
% PURPOSE:  Calculate the left-hand side of the GIM equations used in the minimum		% 
% 	    residual approach.									%
% S. Miller, Oct. 9, 2018									%
%												%
%-----------------------------------------------------------------------------------------------%


function [ A ] = Ax(R, X, HX, D, E, Hpath, weights);

	% FUNCTION INPUTS:
	% R:		Model-data mismatch covariance matrix.
	% X:		Matrix of auxiliary variables.
	% HX:		Auxiliary variables passed through the forward model.
	% D:		Matrix describing temporal covariances in the fluxes.
	% E:		Matrix describing spatial covariances in the fluxes.
	% Hpath:	Path to the H matrix strips.	
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
        m1  = size(E,1);
	ntimes = size(D,1);	

        xi   = weights(1:(size(weights,1)-p));
        betas = weights((size(weights,1)-p+1):size(weights,1));


%----------------------------%
% Calculate the first line   %
%----------------------------%

        %! Note: Edit this section to match the actual format of H in your problem.

        % Calculate H * xi
        % disp('Calculate H^T * xi');
        % tic;
        Hxi = [];
        for j = 1:size(D,1);
        H = load(strcat(Hpath,'H_',num2str(j),'.mat'));
        H = H.H;
	Hxi = [Hxi; H'*xi];
        end;
        % disp(toc);

        % Calculate Q * H^T * xi
        % NOTE: This step is the slowest step for a year-long inverse model.
	% disp('Calculate Q * H^T * xi');
        % tic;
        QHxi = [];

        for j = 1:ntimes;
        A1 = zeros(m1,1);
                for i = 1:ntimes;
                sel = (m1.*(i-1)+1):(i.*m1);
                A1 = A1 + Hxi(sel,1) .* D(j,i);
                end; % End of i loop
        temp = E * A1;
        QHxi = [QHxi; temp];
        end; % End of j loop
        A1 = []; temp =[];
        % disp(toc);

	% Calculate H * Q * H' * xi
        % disp('Calculate HQHxi');
        % tic;
        HQHxi = zeros(n,1);
        for j = 1:ntimes;
        H = load(strcat(Hpath,'H_',num2str(j),'.mat'));
        H = H.H;
	sel   = (m1.*(j-1)+1):(j.*m1);
        HQHxi = HQHxi + H*QHxi(sel,:);
        H = [];
        end;
        % disp(toc);

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

