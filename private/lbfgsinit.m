function H = lbfgsinit(n,k,dscale)
%lbfgsinit  initializes an L-BFGS matrix.
%
%   H = lbfgsinit(n,k,dscale)
%
%   See also lbfgsadd, lbfgsdel, lbfgshprod, lbfgsupdate.

% $Id: lbfgsinit.m 385 2008-09-02 22:22:41Z mpf $

if nargin < 3 || isempty(dscale)
   dscale = 1;
end

H.jNew = 1;
H.jOld = 1;
H.jMax = k;
H.S    =       zeros(n,k);
H.Y    =       zeros(n,k);
H.r    =       zeros(1,k);
H.H0   = dscale*ones(n,1);
