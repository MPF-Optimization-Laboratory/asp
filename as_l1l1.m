function [x,r,inform] = as_l1l1(A, b, lambda, opts, inform)
%AS_L1L1  Solve the one-norm regularized one-norm misfit problem.
%
%   [X,R,INFORM] = AS_L1L1(A,B,LAM) solves the problem
%
%   (L1-L1) minimize_x  LAM ||x||_1  +  ||Ax-B||_1.
%
%   Setting LAM = 0 or LAM = [] (or excluding it) yields a basis
%   pursuit solution.
%
%   AS_L1L1(A,B,LAM,OPTS) specifies options that can be set using
%   AS_SETPARMS.
%
%   AS_L1L1(A,B,LAM,OPTS,INFORM) uses information stored in INFORM
%   (from a previous call to AS_L1L1) to warm-start the
%   algorithm. Note that the previous call to AS_L1L1 must have been
%   to a problem with the same A and B.
%
%   Inputs
%   A       is an m-by-n matrix, explicit or an operator (required).
%   B       is an m-vector (required).
%   LAM     is a nonnegative scalar.
%   OPTS    is an options structure created using AS_SETPARMS.
%   INFORM  is an information structure from a previous call to AS_L1L1.
%
%   Outputs
%   X       is the final estimate of x     (typically including many 0s).
%   R       is the final estimate of b-A*x (typically including many 0s).
%   INFORM  is a structure containing statistics on the solution process.
%
%   Example
%   m = 600; n = 2560; k = 20;    % No. of rows, columns, and nonzeros
%   p = randperm(n); p = p(1:k);  % Position of nonzeros in x
%   x = zeros(n,1);               % Generate sparse solution
%   x(p) = randn(k,1);
%   A = randn(m,n);               % Gaussian m-by-n ensemble
%   b = A*x;                      % Compute the RHS vector
%   [x,r,inform] = as_l1l1(A,b);  % Solve the basis pursuit problem
%
%   See also AS_BPDN, AS_TOPY, AS_SETPARMS, BPDUAL, BPPRIMAL.
%
% ASP Toolbox
% Copyright 2008, Michael P. Friedlander and Michael A. Saunders
% http://www.cs.ubc.ca/labs/scl/asp

% 04 Oct 2008: First version.
%
% Our approach is to recast (L1-L1) as
%
% minimize  ||x||_1 + ||y||_1
% subj to    Ax  +  lam y = b,
%
% and use BPdual to solve its dual:
%
% maximize  b'z - 1/2 sigma ||z||^2
% subj to   ||A'z||_inf <= 1      :     x     (n)
%           ||  z||_inf <= 1/lam  : lam*y = r (m)
%
% We'll actually set sigma = 0 and let BPdual choose a tiny parameter
% to keep z regularized.

% $Id: as_l1l1.m 404 2008-10-05 17:16:48Z mpf $

% Check arguments
if nargin <  2, error('At least 2 arguments needed'); end
if nargin <  3 || isempty(lambda), lambda =  0; end
if nargin <  4, opts = as_setparms; end
if nargin <  5 || isempty(inform)
    [active,state,w,S,R] = deal([]);
else
    active = inform.active;
    state  = inform.state;
    w      = inform.w;
    S      = inform.S;
    R      = inform.R;
end
[m,n] = size(A);

% Setup data for dual problem.
sigma = 0;
AI    = @(x,mode)AIprod(x,mode,A); % AI = [A I]
AIop  = operator(AI,m,m+n);
bl    = [ ones(n,1)                % Bnds on ||A'y||_inf
          repmat(1/lambda,m,1) ];  % Bnds on ||  y||_inf. OK if infinite.

% Fire up BPdual 
[active,state,xx,w,S,R,inform] = ...
    BPdual(AIop,b,-bl,bl,sigma,active,state,w,S,R,opts);

% Extract dual variables, and the residual.
x  = zeros(n,1);
r  = zeros(m,1);
ix = find(active <= n);
ir = find(active >  n);
x(active(ix)  ) = xx(ix);
r(active(ir)-n) = xx(ir);

% If the user wants inform, add warm-start data
inform.lam    = lambda;
inform.w      = w;
inform.active = active;
inform.state  = state;
inform.S      = S;
inform.R      = R;

end % function as_l1l1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRIVATE FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = AIprod(x,mode,A)
[m,n] = size(A);
if mode == 1
    z = A*x(1:n) + x(n+1:n+m);
else
    z(n+1:n+m) = x;  % Simultaneously allocates and assigns
    z(1:n) = A'*x;
end
end % function AIprod
