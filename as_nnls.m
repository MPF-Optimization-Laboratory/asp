function [x,inform] = as_nnls(A, b, c, opts, inform)
%AS_NNLS  Solve the nonnegative least-squares problem.
%
%   [X,INFORM] = AS_NNLS(A,B) solves the problem
%
%     minimize_x  1/2 ||Ax-B||_2^2        subject to  x >= 0.
%
%   AS_NNLS(A,B,C) additionally includes a linear term in the objective:
%
%     minimize_x  1/2 ||Ax-B||_2^2 + C'x  subject to  x >= 0.
%
%   AS_NNLS(A,B,C,OPTS) specifies options that can be set using
%   AS_SETPARMS.
%
%   AS_NNLS(A,B,C,OPTS,INFORM) uses information stored in INFORM
%   (from a previous call to AS_NNLS) to warm-start the
%   algorithm. Note that the previous call to AS_NNLS must have been
%   to a problem with the same A and C.
%
%   In all cases, the INFORM output argument is optional, and contains
%   statistics on the solution process.
%
%   Inputs
%   A       is an m-by-n matrix, explicit or an operator.
%   B       is an m-vector.
%   C       is an n-vector.
%   OPTS    is an options structure created using AS_SETPARMS.
%   INFORM  is an information structure from a previous call to AS_NNLS.
%
%   Example
%   m = 600; n = 2560; k = 20;    % No. of rows, columns, and nonzeros
%   p = randperm(n); p = p(1:k);  % Position of nonzeros in x
%   x = zeros(n,1);               % Generate sparse nonnegative solution
%   x(p) = max(randn(k,1),0);
%   A = randn(m,n);               % Gaussian m-by-n ensemble
%   b = A*x;                      % Compute the RHS vector
%   x = as_nnls(A,b);             % Solve the basis pursuit problem
%
%   See also BPDUAL, AS_SETPARMS.
%
% BPdual Toolbox
% Copyright 2010, Michael P. Friedlander and Michael A. Saunders
% http://www.cs.ubc.ca/labs/scl/bpdual

% 12 Dec 2012: NNLS needs lam = 1 (not 0!).
% $Id: as_nnls.m 524 2010-04-27 22:26:21Z mpf $

% Check arguments
if nargin <  2, error('At least 2 arguments needed'); end
if nargin <  3 || isempty(c), c =  0; end
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

lam = 1;

% Fire up BPdual
[active,state,xx,y,S,R,inform] = ...
    BPdual(A,b,-inf,c,lam,active,state,y,S,R,opts);

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

end % function as_nnls
