function [x,inform] = as_ktopy(A, b, k, opts)
%AS_KTOPY  Apply homotopy to obtain a k-sparse solution.
%
%   [X,INFORM] = AS_KTOPY(A,B,K) applies homotopy to the problem
%
%      (BPDN) minimize_x  1/2 ||Ax - B||_2^2 + LAM ||x||_1
%
%   with a decreasing sequence of LAM. It exits when the iterate is
%   k-sparse.
%
%   AS_KTOPY(A,B,K,OPTS) specifies options that can be set using
%   AS_SETPARMS.
%
%   In all cases, the INFORM output argument is optional, and contains
%   statistics on the solution process.
%
%   Inputs
%   A       is an m-by-n matrix, explicit or an operator.
%   B       is an m-vector.
%   K       is a nonnegative integer.
%   OPTS    is an options structure created using AS_SETPARMS.
%
%   Example
%   m = 600; n = 2560; k = 20;    % No. of rows, columns, and nonzeros
%   p = randperm(n); p = p(1:k);  % Position of nonzeros in x
%   x = zeros(n,1);               % Generate sparse solution
%   x(p) = randn(k,1);
%   A = randn(m,n);               % Gaussian m-by-n ensemble
%   b = A*x;                      % Compute the RHS vector
%   x = as_ktopy(A,b,10);         % Find a 10-sparse estimate of the soln
%
%   See also AS_TOPY, AS_BPDN, AS_SETPARMS, BPDUAL.
%
% BPdual Toolbox
% Copyright 2008, Michael P. Friedlander and Michael A. Saunders
% http://www.cs.ubc.ca/labs/scl/bpdual

% Check arguments
if nargin <  3, error('At least 3 arguments needed'); end
if nargin <  3 || isempty(k), k = inf; end
if nargin <  4, opts = as_setparms; end
    
% Add BPdual options
opts.homotopy = true;

% Install callback.
opts.callback = @(itn,x,active,var,lambda,r,y,pObj,dObj)...
    callback(itn,x,active,var,lambda,r,y,pObj,dObj,k);

% Fire up BPdual 
[active,state,~,y,S,R] = ...
    BPdual(A,b,-1,1,0,[],[],[],[],[],opts);

% Find an LS solution using this active set.
xx = S\b;

% BPdual's solution x is short. Make it full length.
n = length(state);
x = zeros(n,1);
x(active) = xx;

% If the user wants inform, add warm-start data
if nargout == 1, return, end
inform.y      = y;
inform.active = active;
inform.state  = state;
inform.S      = S;
inform.R      = R;

end % function as_ktopy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Private functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function exit = callback(~,~,active,~,~,~,~,~,~,k)
%CALLBACK  Request exit when active set has k or more elements.
exit = length(active) >= k;
end