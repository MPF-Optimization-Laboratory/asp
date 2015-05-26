function y = mtimes(A,B)

% Note that either A or B is an operator since this function gets
% called for both M*C, C*M, where C is the class and M is a matrix or
% vector. This gives the following options, with S for scalar and C
% for any instance of a classOp
%
% 1) M*C, to be implemented as (C'*M')'
% 2) C*M
% 3) s*C
% 4) C*s
% 5) C*C, either of which can be a foreign class

% ASP Toolbox
% Copyright 2008, Michael P. Friedlander and Michael A. Saunders
% http://www.cs.ubc.ca/labs/scl/asp
%
% $Id: mtimes.m 455 2009-05-11 22:05:16Z mpf $
% $Id: mtimes.m 455 2009-05-11 22:05:16Z mpf $

if isnumeric(A)
    % Mode 1 or 3, A*B, with A matrix or scalar
    if isscalar(A) && (size(B,1) ~= 1)
        error('scalar * operator not supported');
    end
    y = (B' * A')';
elseif isnumeric(B)
    % Mode 2 or 4, with A operator and B matrix or scalar
    if isscalar(B) && (size(A,2) ~= 1)
        error('operator * scalar not supported');
    end
    [m,n] = size(A);
    [p,q] = size(B);
    
    % Raise an error when the matrices do not commute
    if p ~= n
        error('Inner matrix dimensions must agree.');
    end
    
    % Pre-allocate y
    y = zeros(m,q);
    
    % Perform operator*vector on each column of B
    for i=1:q
        if A.adjoint == 0
            if isnumeric(A.op)
                y(:,i) = A.op*B(:,i);
            else
                y(:,i) = A.op(B(:,i),1);
            end
        else
            if isnumeric(A.op)
                y(:,i) = A.op'*B(:,i);
            else
                y(:,i) = A.op(B(:,i),2);
            end
        end
    end
else
    % Mode 5. Not supported.
    error('Operation not supported.');
end
