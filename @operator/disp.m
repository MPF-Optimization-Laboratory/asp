function disp(A,name)
%DISP  Display an operator.
%
%   DISP(A) displays an operator, excluding its name.
%
%   DISP(A,NAME) displays an operator along with its name.

% ASP Toolbox
% Copyright 2008, Michael P. Friedlander and Michael A. Saunders
% http://www.cs.ubc.ca/labs/scl/asp
%
% $Id: disp.m 455 2009-05-11 22:05:16Z mpf $

if nargin < 2 || isempty(name)
   name = 'ans';
end

[m,n] = size(A);

fprintf('%s is an operator of size %i x %i\n',name,m,n);
