function [x,r,inform] = as_wls(A, b, d, lambda, opts, inform)
%AS_GL1  Solve the "weighted" regularized LS problem.
%
%   [x,R,INFORM] = AS_WLS(A,b,d,lambda) solves the problem
%
%   (WLS) minimize_x  lambda ||x||_1  +  1/2||D^(-1)(Ax-b)||_2^2,
%
%   where D = diag(d) and d is a vector of nonnegative weights. Setting LAM
%   = 0 or LAM = [] (or excluding it) yields a solution where Ax=b.
%
%   AS_WLS(A,b,d,LAM,OPTS) specifies options that can be set using
%   AS_SETPARMS.
%
%   AS_WLS(A,b,d,lambda,OPTS,INFORM) uses information stored in INFORM
%   (from a previous call to AS_WLS) to warm-start the algorithm. Note that
%   the previous call to AS_WLS must have been to a problem with the same
%   input parameters A,b,d.
%
%   Inputs
%   A       is an m-by-n matrix, explicit or an operator (required).
%   B       is an r-by-s matrix, explicit or an operator (required).
%   b       is an m-vector (required).
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
%   See also AS_BPDN, AS_TOPY, AS_SETPARMS, BPDUAL.
%
% ASP Toolbox
% Copyright 2008, Michael P. Friedlander and Michael A. Saunders
% http://www.cs.ubc.ca/labs/scl/asp

% 27 Jul 2011: First version.

% Check arguments
if nargin <  3, error('At least 2 arguments needed'); end
if nargin <  4 || isempty(lambda), lambda =  0; end
if nargin <  5, opts = as_setparms; end
if nargin <  6 || isempty(inform)
    [active,state,w,S,R] = deal([]);
else
    active = inform.active;
    state  = inform.state;
    w      = inform.w;
    S      = inform.S;
    R      = inform.R;
end
[mA,nA] = size(A);
[mB,nB] = size(B);

% Setup data for dual problem.
sigma = 0;
DABI   = @(x,mode)DABIprod(x, mode, A, B, sigma); % ABI = inv(D)*[A 0; B I]
DABIop = operator(DABIProd, mA+mB, nA+nB);
bnds  = [ repmat( 0, nA, 1)         % x(nA)  gets weight 0
          repmat( 1, mB, 1)         % s(mB)  gets weight 1
        ];

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
function z = DABIprod(x, mode, A, B, sigma)
[m,n] = size(A);
if mode == 1
    z = A*x(1:n) + x(n+1:n+m);
else
    z(n+1:n+m) = x;  % Simultaneously allocates and assigns
    z(1:n) = A'*x;
end
end % function AIprod
