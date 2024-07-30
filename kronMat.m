classdef kronMat
  %kronMat class
  %  A kronMat object is used to represent a matrix A,
  %       where A is a Kronecker product: A = B kron C
  %  Efficient matrix-vector and matrix-transpose-vector
  %       operations are implemented.
  %
  %  The kronMat class has inputs:
  %     B  - First matrix in Kronecker product
  %     C  - Second matrix in Kronecker product
  %   and is based on a structure with three fields:
  %     B
  %     C
  %     transpose - indicates if the matrix has been transposed
  %
  %  Calling Syntax:
  %
  %    A = kronMat    (returns object with empty fields)
  %    A = kronMat(kronMat) (returns kronMat)
  %    A = kronMat(B, C)
  %
  % J. Chung, 1/31/2017
  % J. Chung, 5/1/2017 Modified to work for sum of Kronecker products
  %                     B and C are cells of matrices
  
  properties
    B
    C
    transpose
  end
  
  methods
    function A = kronMat(varargin) % Constructor
      switch nargin
        case 0
          A.transpose = false;
          A.B = [];
          A.C = [];
        case 1
          if isa(varargin{1}, 'kronMat')
            A = varargin{1};
          else
            error('Incorrect input arguments')
          end
        otherwise
          A.transpose = false;
          if nargin == 2
            A.B = varargin{1};
            A.C = varargin{2};
          else
            error('Incorrect number of input arguments')
          end
      end
    end
    
    function A = ctranspose(A) % Overload transpose
      A.transpose = not(A.transpose); % switches booleen transpose flag
    end
    
    function y = mtimes(arg1, arg2) % Overload matrix vector multiplicaiton
      %   Implement A*s and A'*s for kronMat object A and "vector" s.
      
      if isa(arg1,'kronMat')
        % check to see of arg2 is a scalar
        if isa(arg2,'double') && length(arg2) == 1
          error('Matrix times a scalar not implemented yet')
        else
          % if numel(arg2) ~= length(arg2) %  Check to see if arg2 is not a vector
          %  warning('Warning:kronMat/mtimes assumes input is a vector')
          % end
         
	  % IF A AND B ARE CELLS 
          if iscell(arg1.B) && numel(arg2) == length(arg2) % If B is a cell and arg2 is a vector
            [m1, m2] = size(arg1.C{1});
            [n1, n2] = size(arg1.B{1});
            
            if arg1.transpose
              % Matrix transpose times vector
              y = zeros(m2,n2);
              for i = 1:length(arg1.B)
                y = y + arg1.C{i}'*reshape(arg2,m1,n1)*arg1.B{i};
              end
              
            else
              % Matrix times vector
              y = zeros(m1,n1);
              for i = 1:length(arg1.B)
                y = y + arg1.C{i}*reshape(arg2,m2,n2)*arg1.B{i}';
              end
            end
            y = y(:);
          
	  % IF ARG2 IS A VECTOR  
          elseif numel(arg2) == length(arg2) % If arg2 is a vector
            [m1, m2] = size(arg1.C);
            [n1, n2] = size(arg1.B);
            
            if arg1.transpose
              % Matrix transpose times vector
              y = arg1.C'*reshape(arg2,m1,n1)*arg1.B;
            else
              % Matrix times vector
              y = arg1.C*reshape(arg2,m2,n2)*arg1.B';
            end
            y = y(:);

	  % IF ARG2 IS A MATRIX (NOT A VECTOR)
          else; % If arg2 is a matrix, not a vector
            [m1, m2] = size(arg1.C);
            [n1, n2] = size(arg1.B);
            p = size(arg2,2); 
            y = zeros(m1.*n1,p);
	    for i = 1:p;
              if arg1.transpose
                % Matrix transpose times vector
                ytemp = arg1.C'*reshape(arg2(:,i),m1,n1)*arg1.B;
              else
                % Matrix times vector
                ytemp = arg1.C*reshape(arg2(:,i),m2,n2)*arg1.B';
              end
	      y(:,i) = ytemp(:);
            end; % End p loop

          end
        end
        
      elseif (( isa(arg1, 'double')) && (length(arg1)==1))
        error('Multiplication is not implemented.')
      else
        error('Multiplication is not implemented.')
      end
    end
    
    function varargout = size(A, dim) % Overload size
      if iscell(A.B)
        d = size(A.B{1}).*size(A.C{1});
      else
        d = size(A.B).*size(A.C);
      end
      if nargin == 2
        d = d(dim);
      end
      
      if nargout == 1 || nargout == 0
        varargout{1} = d;
      else
        for i = 1:length(d)
          varargout{i} = d(i);
        end
      end
      
    end
    
  end % methods
  
  
  
end % classdef
