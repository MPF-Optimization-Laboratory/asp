function result = ctranspose(A)
%'  Complex conjugate tranpose.
%   A' is the complex conjugate transpose of A.
%
%   CTRANSPOSE(A) is called for the syntax A' when A is an operator.

% ASP Toolbox
% Copyright 2008, Michael P. Friedlander and Michael A. Saunders
% http://www.cs.ubc.ca/labs/scl/asp
%
% $Id: ctranspose.m 455 2009-05-11 22:05:16Z mpf $

A.adjoint = xor(A.adjoint,1);

result = A;
