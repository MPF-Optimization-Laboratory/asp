function x = linsolve(A,b,opts)
%LINSOLVE  returns the result of  A\b.
%
%   This routine is simply a wrapper to Octave's back-slash operator.
%   Its purpose is to make ASP compatible with Octave.
%
% BPdual Toolbox
% Copyright 2008, Michael P. Friedlander and Michael A. Saunders
% http://www.cs.ubc.ca/labs/scl/bpdual

% $Id: linsolve.m 389 2008-09-03 06:03:53Z mpf $

if isfield(opts,'TRANSA') && opts.TRANSA
    x = A'\b;
else
    x = A\b;
end
