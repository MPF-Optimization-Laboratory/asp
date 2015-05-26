function [t, k, stat, ix, ex] = l2l1line(r,dr,x,dx,alpha,beta,gamma)

% l2l1line solves the univariate problem
%
% min_t  q(t) =   1/2 alpha||r + t dr||^2
%               + 1/2 beta ||x + t dx||^2
%               +     gamma||x + t dx||_1
%

% 30 Oct 2007: First version of l2l1line.
%              q(t) = 1/2 ||r + t dr||^2 + lam ||x||_1
% 09 Aug 2008: Modified to allow for a 2-norm term in x.
%
% $Id: BPdual.m 202 2007-09-14 20:43:21Z saunders $
%----------------------------------------------------------------------

debug    = false;                       % Set to 'true' to print log
n        = length(x);
h        = alpha*(dr'*dr) + beta*(dx'*dx); % Directional Hessian (constant)
STEP_LHS = 1;
STEP_RHS = 2;
STEP_MID = 3;
STEP_INF = 4;
BigNum   = 1e20;
ex       = [];
ix       = [];

% pivTol is used to determine if a step is nontrivial.
% Scale pivTol if norm(dx) is large.
pivTol = 1e-10*max(1, norm(dx,inf));

% ----------------------------------------------------------------------
% Compute the list of sorted breakpoints.
% ----------------------------------------------------------------------
inc = dx >  pivTol  &  x <= 0;
dec = dx < -pivTol  &  x >= 0;
mov = inc | dec;

aBrk      = inf(n+1,1);
aBrk(mov) = - x(mov) ./ dx(mov);

[aBrk, iBrk] = sort(aBrk, 'ascend');

% ----------------------------------------------------------------------
% Compute directional derivative at t = 0:
%
%   q'(0) = r'*dr + lam * e' * dx  where  e = sgn(x).
%
% Special care must be taken to define sgn(x) for variables that are on or
% very near their breakpoints.  For these variables, sgn(x) is determined
% by the sign of the corresponding search direction.
% ----------------------------------------------------------------------
e = sign(x);
k = 1;
while aBrk(k) < pivTol
    ix    = iBrk(k);
    e(ix) = sign( dx(ix) );
    k     = k + 1;
    if k > n, break; end
end
g = alpha*(r'*dr) + beta*(x'*dx) + gamma*(e'*dx); % Update directional deriv.

% ----------------------------------------------------------------------
% Print to the log.
% ----------------------------------------------------------------------
if debug
   nBrk = length(find(mov));
   I = iBrk(k:nBrk);
   fprintf('\n\n=========== Begin linesearch ================\n');
   fprintf('\n %5s %10s %10s %4s %10s\n',...
       'ix','x','dx','e','break');
   fprintf(' %5i %10.3e %10.3e %4i %10.3e\n',...
       [I,x(I),dx(I),e(I),aBrk(k:nBrk)]');
   fprintf('\n %4s %4s %10s %10s %10s %5s %10s %10s\n',...
       'k','ix','step','t','dx','e','g','q')
end

% ----------------------------------------------------------------------
% Explore each breakpoint.  There are at most n of them.
% ----------------------------------------------------------------------
t = 0; rhs = 0;
while true
    
    if g >= 0
       % Minimum is on LHS. Done! (Captures the dr=dx=0 case.)
       stat = STEP_LHS;
       step = 0;
       
    else
       % Minimizer must be to the right.
       % Grab the boundaries of the current search interval.
       lhs     = rhs;
       rhs     = aBrk(k);       % Next breakpoint...
       ix      = iBrk(k);       % and the corresponding variable.
       maxStep = rhs - lhs;     % Maximum step to that breakpoint.
    
       if h < eps
          % If we got this far, g negative and q is nearly linear.
          % Minimizer is either at infinity or on the RHS.
          if maxStep >= BigNum
             stat = STEP_INF;
          else
             stat = STEP_RHS;
             step = maxStep;
          end
        
       else
          % Minimum is on RHS or inside the interval. Check.
          step = -g / h;
          if step >= BigNum  &&  step < maxStep
             % Step is unbounded.  Exit with unbounded flag.
             stat = STEP_INF;
          elseif step >= maxStep
             % Minimum is on RHS.  Truncate the step and continue.
             stat = STEP_RHS;
             step = maxStep;
             t    = rhs;
          else
             % Minimum found in the interval's interior.  Exit.
             stat = STEP_MID;
             t    = t + step;
          end 
       end
    end
    
    if debug
       q = .5*alpha*norm(r+t*dr)^2 + .5*beta*norm(x+t*dx)^2 + gamma*norm(x+t*dx,1);
       fprintf(' %4i %4i %10.3e %10.3e %10.3e %5i %10.3e %14.7e',...
           k,ix,step,t,dx(ix),e(ix),g,q)
       if     stat == STEP_RHS, fprintf(' %s\n','RHS');
       elseif stat == STEP_MID, fprintf(' %s\n','MID');
       elseif stat == STEP_INF, fprintf(' %s\n','INF');
       elseif stat == STEP_LHS, fprintf(' %s\n','LHS');
       end
    end
    
    % -----------------------------------------------------------------
    % Check for exit.
    % -----------------------------------------------------------------
    if stat == STEP_RHS
       % We hit a breakpoint.  Keep going!
    else
       break
    end
    
    % -----------------------------------------------------------------
    % We'll now explore the next interval. Update the gradient.
    % -----------------------------------------------------------------
    ex = e(ix);
    e(ix) = -ex;
    g = g + 2 * gamma * e(ix) * dx(ix);  % passed through a kink
    g = g + step * h;                    % account for the 2nd-order term
    k = k + 1;

end  % while true

% k is currently the RHS breakpoint. We'll return the LHS breakpoint.
k = k - 1;

% We should never exit with a STEP_RHS status.
if stat == STEP_RHS
   error('Incorrect stat in linesearch')
end

if debug,
    fprintf('\n=========== end linesearch ================\n');
end

end % function
