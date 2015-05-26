function [x,inform] = as_topy(A, b, lam, opts, inform)
%AS_TOPY  Apply the Homotopy method to the basis pursuit problem.
%
%   [X,INFORM] = AS_TOPY(A,B) solves the basis pursuit problem
%
%      (BP)   minimize_x  ||x||_1  subject to  Ax = B.
%
%   AS_TOPY(A,B,LAM) solves the basis pursuit denoise problem
%
%      (BPDN) minimize_x  1/2 ||Ax - B||_2^2 + LAM ||x||_1.
%
%   Setting LAM = 0 or LAM = [] will yield exactly the same solution
%   as the call AS_TOPY(A,B).
%
%   AS_TOPY(A,B,LAM,OPTS) specifies options that can be set using
%   AS_SETPARMS.
%
%   AS_TOPY(A,B,LAM,OPTS,INFORM) uses information stored in INFORM
%   (from a previous call to AS_TOPY) to warm-start the
%   algorithm. Note that the previous call to AS_TOPY must have been
%   to a problem with the same A and B.
%
%   In all cases, the INFORM output argument is optional, and contains
%   statistics on the solution process.
%
%   Inputs
%   A       is an m-by-n matrix, explicit or an operator.
%   B       is an m-vector.
%   LAM     is a nonnegative scalar.
%   OPTS    is an options structure created using AS_SETPARMS.
%   INFORM  is an information structure from a previous call to AS_TOPY.
%
%   Example
%   m = 600; n = 2560; k = 20;    % No. of rows, columns, and nonzeros
%   p = randperm(n); p = p(1:k);  % Position of nonzeros in x
%   x = zeros(n,1);               % Generate sparse solution
%   x(p) = randn(k,1);
%   A = randn(m,n);               % Gaussian m-by-n ensemble
%   b = A*x;                      % Compute the RHS vector
%   x = as_topy(A,b);             % Solve the basis pursuit problem
%
%   See also AS_BPDN, AS_SETPARMS, BPDUAL.
%
% BPdual Toolbox
% Copyright 2008, Michael P. Friedlander and Michael A. Saunders
% http://www.cs.ubc.ca/labs/scl/bpdual

% Check arguments
if nargin <  2, error('At least 2 arguments needed'); end
if nargin <  3 || isempty(lam), lam  =  0; end
if nargin <  4, opts = as_setparms; end
if nargin <  5 || isempty(inform)
    [active,state,y,S,R] = deal([]);
else
    active = inform.active;
    state  = inform.state;
    y      = inform.y;
    S      = inform.S;
    R      = inform.R;
end
if ~isscalar(lam)
    error('LAM must be a scalar.');
end

% Add BPdual options
opts.homotopy = true;

% Fire up BPdual 
[active,state,xx,y,S,R,inform] = ...
    BPdual(A,b,-1,1,lam,active,state,y,S,R,opts);

% BPdual's solution x is short. Make it full length.
n = length(state);
x = zeros(n,1);
x(active) = xx;

% If the user wants inform, add warm-start data
if nargout == 1, return, end
inform.lam    = lam;
inform.y      = y;
inform.active = active;
inform.state  = state;
inform.S      = S;
inform.R      = R;

end % function as_topy
