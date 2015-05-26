function y = aprodind(mode,x,A,m,n,idx)
%aprodind  matrix-vector multiplication with a submatrix of Aprod
%
%   Y = APRODIND(MODE,X,A,M,N,IDX) returns the product
%
%   Y = A(:,IDX) *X   if  MODE == 1,
%   Y = A(:,IDX)'*X   if  MODE == 2.
%
%   where A is an M-by-N matrix or operator.
%
%BPdual Toolbox
%Copyright 2008, Michael P. Friedlander and Michael A. Saunders
%http://www.cs.ubc.ca/labs/scl/bpdual

%$Id: as_rwbp.m 343 2008-08-01 21:34:53Z mpf $

if mode == 1
   z = zeros(n,1);
   z(idx) = x;
   y = A*z;
else
   z = A'*x;
   y = z(idx);
end
