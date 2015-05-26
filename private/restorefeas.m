function y = restorefeas(y,active,state,S,R,bl,bu)
%RESTOREFEAS  Restore feasibility of y with respect to active constraints.
%
% Find a least-norm change dy to y so that y+dy is feasible for this active
% set.  Thus, solve the least squares problem
%
%   minimize  ||dy||
%   subj to   S'(y + dy) = c
%
% where  c(i) = bl(i) if state(i) == -1
%             = bu(i) if state(i) == +1.
%
% This problem has optimality conditions
% [ -I   S ] [dy] = [0]
% [  S'    ] [ x]   [c - S'y],
% which we solve in two steps:
% 1. R'R x = c - S'y  (where R'R = S'S)
% 2. dy = Sx.
%
% NOTE: The final value of z = A'y may violate inactive constraints.

% 30 Jul 2008: First version.
% 16 Aug 2009: Guard against infinite bounds.  Inf .* false give Nan!

lbnd = state(active) == -1 & bl(active) > -inf;
ubnd = state(active) == +1 & bu(active) < +inf;
c = bl(active).*lbnd + bu(active).*ubnd;

c  = c - S'*y;
x  = R' \ c;
x  = R  \ x;
dy = S*x;
y  = y + dy;
  
end % function restorefeas
