function [R,noup] = bfgsupdate(R, step, p, g, g2, v)
%bfgsupdate  returns a BFGS update to the triangular factor R.
%
%   [R,noup] = bfgsupdate(R, step, p, g, g2, v) returns a BFGS update
%   to the triangle R.  The flag noup indicates if the updated was
%   skipped.
%
%   On entry, v solves R'v = g.

% $Id$

gtp  = g' *p;
gtp2 = g2'*p;
noup = gtp2 <= 0.91*gtp; 
if noup
   return
end

delta1 = 1 / sqrt( abs( gtp ) );
delta2 = 1 / sqrt( step*(gtp2 - gtp) );

w  = - delta1*v;
z  =   delta2*(g2-g) + delta1*g;
R  =   rank1triu(R,w,z);
