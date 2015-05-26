function A = operator(op,m,n)
%OPERATOR  Create a linear operator.
%
%   A = OPERATOR(OP,M,N) creates a M-times-N matrix-like linear operator
%   from the linear function
%
%   y = OP(X,MODE).
%
%   The resulting object A can be used like a regular Matlab matrix, i.e.,
%
%        y = A *x  is equivalent to  y = OP(x,1)
%   and  z = A'*y  is equivalent to  z = OP(y,2).
%
%   See also @OPERATOR/MTIMES, @OPERATOR/SIZE, @OPERATOR/CTRANSPOSE.
%
% ASP Toolbox
% Copyright 2008, Michael P. Friedlander and Michael A. Saunders
% http://www.cs.ubc.ca/labs/scl/asp

% $Id: operator.m 455 2009-05-11 22:05:16Z mpf $

% Empty constructor
if nargin == 0
    A.adjoint = 0;
    A.op      = [];
    A.size    = [0 0];
    A         = class(A,'operator');
    return;
end

% Conversion constructors.
if isa(op,'operator')
    A = op;
   
elseif isa(op,'function_handle')
    if nargin < 3
        error('Missing dimension arguments.');
    end
    A.adjoint = 0;
    A.op      = op;
    A.size    = [m n];
    A         = class(A,'operator');
    
elseif isnumeric(op)
    A.adjoint = 0;
    A.op      = op;
    A.size    = size(op);
    A         = class(A,'operator');
    
else
    error('Unsupported use of function OPERATOR.');

end
