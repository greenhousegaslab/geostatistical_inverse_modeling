%-----------------------------------------------------------------------------------------------------------------------------------------------%
% FUNCTION: bic.m																%
% PURPOSE:  Compute the BIC scores for all possible candidate models. This script is specifically designed for analyzing CARVE methane data.	%
%	Do NOT use branch and bound														%
% S. Miller, July 15, 2015															%
% 	Updated Aug. 28, 2025
%																		%
%-----------------------------------------------------------------------------------------------------------------------------------------------%

	
	%-------------%
	% Notes:      %
	%-------------%
		
	
function [ selmodel beta1 ] = bic(Z, A, B, C, p, detpsi);

	% FUNCTION INPUTS:
	% 	p is the number of covariates (excluding the intercept). The intercept should be the first columns/rows in B and C

	

%------------------------------------%
% 1. Create initial model set M_p    %
%------------------------------------%

	disp('Evaluate the initial model set');

        % Set the number of intercept terms
        p0 = size(C,1) - p;
	n  = length(Z);

        % Create the final output cell
        output = cell(0,4);
	
	
%--------------------------------------%
% Loop over models of different sizes  %
%--------------------------------------%

	for nvar = 1:p;

	disp('Considering candidate models with this number of variables:');
	disp(num2str(nvar));
	
	% Create the list of possible models
	varlist = combnk(1:p,nvar);
	
	% Create an index for all models under consideration
	modelindex = 1:size(varlist,1);
	disp('Number of combinations');
	disp(length(modelindex));
	
	% Put the combination into a cell structure
        combos         = cell(size(varlist,1),4);
        % Column 1: covariates
        % Column 2: RSS
        % Column 3: BIC score
	% Column 4: Drift coefficients
        for j = modelindex;
        combos{j,1} = varlist(j,:);
        end;


%---------------------------------------%
% 2. Evaluate the initial model set     %
%---------------------------------------%

        % Create a temporary output
        tempout = cell(size(varlist,1),3);

        % Loop over each model
        for j = modelindex;

        % Select out the appropriate elements of B and C
        % Account for the intercept term
        sel = [ 1:p0 (combos{j,1}+p0)];
        Btemp = B(:,sel);
        Ctemp = C(sel,sel);

	beta1 = Ctemp \ Btemp';
	RSS   = A - Btemp * (Ctemp \ Btemp'); 

        % Calculate the BIC score
        % Note that BIC calculation typically includes log(det(Psi)). 
        % For now, we've left out that component because it is constant for our purposes 
        % (b/c we're not varying the covariance matrix parameters).
        BIC = RSS + length(sel) .* log(n);

        % Put the output into the temporary matrix
        tempcell = cell(1,3);
        tempcell{1} = RSS;
        tempcell{2} = BIC;
	tempcell{3} = beta1;
        tempout(j,:) = tempcell;

        end; % End of the j loop

        % Merge the temporary output cell with the full output cell
        combos(:,2:4) = tempout;
	
	% Append 'combos' to the final output array
	output = [output; combos];

	end; % End of nvar loop


%----------------------------------%
% 8. Gather up the final results   %
%----------------------------------%

        % Sort the outputs by BIC score
        output = sortrows(output,3);
	
	% Find the first model that doesn't have negative coefficients on the auxiliary variables
	
	selmodel   = 0;
	selcounter = 1;
	while selmodel == 0;
	coeff = output{selcounter,4};
	mincoeff = min(coeff(2:length(coeff)));
	if mincoeff > 0; selmodel = selcounter; end;
	selcounter = selcounter + 1;
	end;

	disp('Number of models with negative coefficients before reaching the selected model');
	disp(num2str(selmodel));

	disp('RSS of the selected model');
	disp(output{selmodel,2});

	disp('BIC score of the selected model');
	disp(output{selmodel,3} + detpsi);

	disp('Drift coefficients of the selected model');
	beta1 = output{selmodel,4};
	disp(beta1);

        % Pull out the index of the covariates
        selmodel = output{selmodel,1};


%-----------------------------------------------------------------------------------------------------------------------------------------------%
% END OF FUNCTION																%
%-----------------------------------------------------------------------------------------------------------------------------------------------%

