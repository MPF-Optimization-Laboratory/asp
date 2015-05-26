function [x, r] = csne(R,A,b)
%CSNE   Solve the corrected semi-normal equations R'Rx=A'b.
%
% x = csne(R,A,b) solves the least-squares problem
%
%   minimize  || Ax - b ||_2
%
% using the corrected semi-normal equation approach described by
% Bjork (1987).
%
% [x,r] = csne(R,A,b) additionally returns the residual r = b - Ax.

% Michael Friedlander (mpf@cs.ubc.ca) and
% Michael Saunders (saunders@stanford.edu)

% linsolve options
optsR.UT       = true;           % for solves with R
optsRT.UT      = true;           % for solves with R'
optsRT.TRANSA  = true;           % ...             R'

q      = A'*b;
x      = linsolve(R,q,optsRT);   % R'x = A'b

bnorm2 = norm(b)^2;
xnorm2 = norm(x)^2;
d2     = bnorm2 - xnorm2;

x      = linsolve(R,x,optsR );   % R x = x

% We hope to save on the iterative refinement step if this tests true.
% See QRaddcol.
if 0%d2 > 0.01*bnorm2
   if nargout > 1
      r = b - A*x;
   end
   return
else
   % Apply one step of iterative refinement.
   r   = b - A*x;
   q   = A'*r;
   dx  = linsolve(R,q ,optsRT);  % R'dx = q
   dx  = linsolve(R,dx,optsR );  % R dx = dx
   x   = x + dx;
end