function D = bfgsdiag(D, step, p, g, g2, dscale)

% bfgsdiag  returns a BFGS update to just the diagonals D
% of the classical BFGS quasi-Newton update to a Hessian:
%   B := B + gg'/g'p + yy'/y's,   where s = step*p,  y = g2-g.
%
% On entry, p satisfies Dp = -g, but we don't use this.

% 20 Aug 2008: First version.
% $Id$

gtp  = g' *p;
gtp2 = g2'*p;

if gtp2 <= 0.91*gtp
  fprintf('\n  W: bfgsdiag update skipped')
  D(:) = dscale;
  return
end

delta1 = 1/gtp;                   % Negative
delta2 = 1/(step*(gtp2 - gtp));   % Positive
D      = D + delta1*g.^2 + delta2*(g2-g).^2;

if any(D<=0)
  fprintf('\n  W: bfgsdiag not posdef')
  D(:) = dscale;
end
