classdef matvecH
    %
    % matvecH class
    %
    % A matvecH object is used to represent a matrix Hall in matrix-vector
    % multiplication where the matrix is stacked by H_i.mat horizontally
    %
    % the matvecH has input(s):
    %   n - number time periods in the inverse model
    %   path - path where H matrices are located
    %
    % and is based on a structure with the following fields:
    %
    % Calling Syntax:
    % P = matvecH(n,path)
    %
    % T.Cho, 10/29/2020
    % Modified by S. Miller, 7/24/2024
    
    properties
        n
        path
        transpose
    end % properties
    
    methods
        
        function P = matvecH(varargin) % constructor
            switch nargin
                case 2
                    P.transpose = false;
                    P.n = varargin{1};
                    P.path = varargin{2};
                otherwise
                    error('Incorrect number of input arguments')
            end % switch
        end % constructor
        
        function P = ctranspose(P) % Overload transpose
            P.transpose = not(P.transpose); % switches boolean transpose flag
        end % transpose
        
        function y = mtimes(A,x)
            load(strcat(A.path,'H_1.mat'));
            [mA,nA] = size(H);
            [mx,nx] = size(x);
            if A.transpose % transpose
                if mx ~= mA
                    error('Invalid size of x') 
                end
                Z = zeros(A.n*nA,nx);
                for i = 1:nx
                   z = zeros(A.n*nA,1); 
                   for j = 1:A.n
                       load(strcat(A.path,'H_',num2str(j),'.mat'));
                       z((j-1)*nA+1:j*nA ,1) = H'*x(:,i);
                   end
                   Z(:,i) = z;
                end     
            else % no transpose
                if mx ~=  A.n*nA
                    error('Invalid size of x') 
                end
                Z = zeros(mA,nx);
                for i = 1:nx
                    X = reshape(x(:,i),[nA A.n]);
                    z = zeros(mA,1);
                    for j = 1:A.n
                        load(strcat(A.path,'H_',num2str(j),'.mat'));
                        z = z+H*X(:,j);
                    end
                    Z(:,i) = z;
                end
            end 
            y = Z;  
        end % mtimes
        
        function varargout = size(A,dim)
            load(strcat(A.path,'H_1.mat'));
            [mA,nA] = size(H);
            d(1) = mA;
            d(2) = nA*A.n;
            if nargout == 1 || nargout == 0
                if nargin >1
                    varargout{1} = d(dim);
                else
                    varargout{1} = d;
                end
            else
                varargout{1} = d(1);
                varargout{2} = d(2);
            end
        end % size
        
        function l = length(A)
            load(strcat(A.path,'H_1.mat'));
            [mA,nA] = size(H);
            if nA*A.n > mA
                l = nA*A.n;
            else
                l = mA; 
            end
        end % length
        
    end % methods
    
end % end classdef