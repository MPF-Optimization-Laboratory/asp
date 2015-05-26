function [x,inform] = as_bpdn(A, b, lam, opts, inform)
%AS_BPDN  Applies the BPdual method to the basis pursuit problem.
%
%   [X,INFORM] = AS_BPDN(A,B) solves the basis pursuit problem
%
%      (BP)   minimize_x  ||x||_1  subject to  Ax = B.
%
%   [X,INFORM] = AS_BPDN(A,B,LAM) solves the basis pursuit denoise problem
%
%      (BPDN) minimize_x  1/2 ||Ax - B||_2^2 + LAM ||x||_1.
%
%   Setting LAM = 0 or LAM = [] will yield exactly the same solution
%   as the call AS_BPDN(A,B).
%
%   AS_BPDN(A,B,LAM,OPTS) specifies options that can be set using
%   AS_SETPARMS.
%
%   AS_BPDN(A,B,LAM,OPTS,INFORM) uses information stored in INFORM
%   (from a previous call to AS_BPDN) to warm-start the
%   algorithm. Note that the previous call to AS_BPDN must have been
%   to a problem with the same A and B.
%
%   In all cases, the INFORM output argument is optional, and contains
%   statistics on the solution process.
%
%   Inputs
%   A       is an m-by-n matrix, explicit or an operator.
%           If A is a function, then it must have the signature
%
%           y = A(x,mode)   if mode == 1 then y = A x  (y is m-by-1);
%                           if mode == 2 then y = A'x  (y is n-by-1).
%
%   B       is an m-vector.
%   LAM     is a nonnegative scalar.
%   OPTS    is an options structure created using AS_SETPARMS.
%   INFORM  is an information structure from a previous call to AS_BPDN.
%
%   Example
%   m = 600; n = 2560; k = 20;    % No. of rows, columns, and nonzeros
%   p = randperm(n); p = p(1:k);  % Position of nonzeros in x
%   x = zeros(n,1);               % Generate sparse solution
%   x(p) = randn(k,1);
%   A = randn(m,n);               % Gaussian m-by-n ensemble
%   b = A*x;                      % Compute the RHS vector
%   [x,inform] = as_bpdn(A,b);    % Solve the basis pursuit problem
%
%   See also AS_TOPY, AS_SETPARMS, BPDUAL.
%
%BPdual Toolbox
%Copyright 2008, Michael P. Friedlander and Michael A. Saunders
%http://www.cs.ubc.ca/labs/scl/bpdual

%$Id: as_bpdn.m 489 2009-08-20 16:43:09Z mpf $

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
if length(lam) > 1, error('LAM (3rd arg) must be a scalar'); end
    
% Fire up BPdual 
[active,state,xx,y,S,R,inform] = ...
    BPdual(A,b,-1,1,lam,active,state,y,S,R,opts);

% BPdual's solution x is short. Make it full length.
x = zeros(size(state));
x(active) = xx;

% If the user wants inform, add warm-start data
inform.lam    = lam;
inform.y      = y;
inform.active = active;
inform.state  = state;
inform.S      = S;
inform.R      = R;
