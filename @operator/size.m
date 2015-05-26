function varargout = size(A,dim)
%SIZE  Dimensions of an operator.
%
%   D = SIZE(A), for an M-by-N an operator A, returns the
%   two-element row vectors D = [M,N].
%
%   [M,N] = SIZE(A) returns M and N as separate arguments.
%
%   M = SIZE(A,DIM) retuns the length of the dimension specified by
%   the scalar DIM.  Note that DIM must be 1 or 2.

% ASP Toolbox
% Copyright 2008, Michael P. Friedlander and Michael A. Saunders
% http://www.cs.ubc.ca/labs/scl/asp
%
% $Id: size.m 455 2009-05-11 22:05:16Z mpf $

sz = A.size;        % Size of the operator.
if A.adjoint
    sz = sz([2 1]);   % We're working with its adjoint. Invert.
end

if nargin > 1
    if nargout > 1
       error('Too many output arguments.');
    end
    
    if ~( dim > 0 && isnumeric(dim) && rem(dim,1) == 0 )
       error(['Dimension argument must be a positive integer' ...
              'scalar within indexing range.']);
    elseif dim > 2
        varargout{1} = 1;
    else
        varargout{1} = sz(dim);
    end
else
    if nargout < 2
        varargout{1} = sz;
    else
        varargout{1} = sz(1);
        varargout{2} = sz(2);
        for i=3:nargout
          varargout{i} = 1;
        end        
    end
end
