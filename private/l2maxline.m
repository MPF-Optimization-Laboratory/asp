function [tStep,state,k,iz,tBreak] = ...
    l2maxline(g,dy,dz,z,bl,bu,lam,rho,state,angles)

% Linesearch on the piecewise-quadratic function
%
% q(t) = t c'dy + 0.5 lam t^2 dy'dy
%        + rho sum[ max(0,  bl - (z + t dz) ) ]
%        + rho sum[ max(0, -bu + (z + t dz) ) ]
%
% The gradient is given by
%
% q'(t) = c'dy + rho e'dz + lam t dy'dy.
%
% where
%
% e(i) = -1  if  state(i) = -2        violated below
%      = +1  if           = +2        violated above
%      =  0  if           = -1,0,+1   active(+/-1) or feasible(0)
%
% and  state  is defined below. 
%
% The first term is independent of t, and so it's only necessary to
% update the 2nd two terms (note that e depends on t).
%
% INPUTS
%
% g: scalar.
%    On entry, it is assumed that
%      g = c'dy + rho e'dz
%    where e is as def'd above and corresponds the value of state on input.
%
% dy, dz: m- and n-vectors.
%    Search direction.
%
% bl, bu: n-vectors.
%    Lower and upper bounds on the (implicit) variable z, which is
%    assumed to be 0, ie, the original bounds on z are really
%    l <= z + t dz <= u, but we recenter, so that
%    bl := l - z <= t dz <= bu := u - z.
%
% lam,rho: scalars.
%
% state: n-vector.
%    Gives the current state of each variable relative to bl,bu.  The
%    value of state(iz) determines if variable iz is currently below,
%    at, or below its bound:
%
%                     ----------bl----------bu----------
%           state =      -2     -1    0     +1    +2
%
%    Assume:
%    1. If state = +2,0,-2, then z is at least fixTol away from its bound.
%    2. Fixed variables (where bu - bl <= fixTol) don't have state = 0.
%
% OUTPUTS
% tStep   is the final steplength
% state   is the final state of each variable after the linesearch.
% k       is the number of breakpoints encountered.  This is also the number
%         of variables that changed state.
% iz      If the final step is constrained, iz is the index of the active
%         breakpoint.  If the final step is unconstrained, iz=0.
% tBreak  = +1 if an upper-bound breakpoint is active
%           -1 if a  lower-bound breakpoint is active
%            0 if the step is unconstrained.
%
% 30 Oct 2007: First version of l2l1line.
% $Id: l2maxline.m 511 2010-04-16 04:02:43Z mpf $
%----------------------------------------------------------------------

debug     = 0;                % 0=no log, 1=itn log, 2=too much log!
backtrack = nargin >= 10;
n         = length(dz);
h         = lam * (dy'*dy);   % 2nd derivative.
STEP_LHS  = 1;
STEP_RHS  = 2;
STEP_MID  = 3;
STEP_INF  = 4;
BigNum    = 1e20;
using_octave = exist('OCTAVE_VERSION','builtin');

% Storage for sequence of encountered breakpoint:
s.viz     = zeros(n,1);  % variable indices
s.state   = state;       % their states
s.k       = 0;           % index of "best" breakpoint
s.aiz     = 0;           % magnitude of "largest" encountered angle
s.tBreak  = nan;         % type of break: -1=lower, +1=upper

% pivTol is used to determine if a step is nontrivial.
% Scale pivTol if norm(dz) is large.
pivTol = 1e-10*max(1, norm(dz,inf));
dzt    = dz .* ( abs(dz) > pivTol );

% Determine if bl(i) =~ bu(i).
fixTol = sqrt(eps);
fixed  = bu - bl < fixTol;

% Misc variables for the log.
logH1 = '\n %4s %10s %10s %10s %10s %10s %10s %5s %5s %5s %5s\n';
logB1 = ' %4i %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %5i %5i %5i %5i\n';
logH2 = '\n %4s %4s %5s %5s %10s %10s %10s %10s %10s %10s %10s %10s\n';
logB2 = ' %4i %4i %5s %5s %10.2e %10.2e %10.2e %+10.2e %+10.2e %+10.2e %+10.2e %+10.2e\n';

% ----------------------------------------------------------------------
% Find the breakpoints.  There are at most 2n breakpoints.
% Lower breakpoints (corresponding to bl) are in 1  :n,
% upper breakpoints (corresponding to bu) are in 1+n:2n.
% ----------------------------------------------------------------------
incl = dz >  pivTol  &  state == -2;
decl = dz < -pivTol  & (state == +2 | state == 0);
movl = find(incl | decl);
sl   = bl - z;

incu = dz >  pivTol  & (state == -2 | state == 0);
decu = dz < -pivTol  &  state == +2;
movu = find(incu | decu);
su   = bu - z;

nBrk = 2*n;           % Potential number of breakpoints
iBrk = (1:nBrk)';
aBrk = inf(nBrk,1);
aBrk(movl  ) = sl(movl) ./ dz(movl);
aBrk(movu+n) = su(movu) ./ dz(movu);

if using_octave
    [aBrk, iBrk] = heapify(aBrk, iBrk);
else
    heapify(aBrk, iBrk);  % !!! This call modifies the RHS args.
end

% ----------------------------------------------------------------------
% Plot the curve if in debug mode.
% ----------------------------------------------------------------------
if debug
    aBrkTmp(iBrk) = aBrk;  % Restore original order.
    [~,iTmp] = sort(aBrkTmp,'ascend');
    iBrkTmp(iTmp) = (1:2*n); clear iTmp;
    fprintf(logH1,'iz','bl','z','dz','bu','breakL','breakU','state','sortL',...
            'sortU','state');
    fprintf(logB1,[(1:n)' bl z dz bu aBrkTmp(1:n)' aBrkTmp(1+n:end)'...
                    state iBrkTmp(1:n)' iBrkTmp(n+1:end)' state]');
end

% ----------------------------------------------------------------------
% Explore each breakpoint.  There are at most n of them.
% ----------------------------------------------------------------------
k = 0; cStep = 0; tStep = 0; tBreak = 0; % The urrent and total steplengths
while nBrk > 0
    
    if g >= 0
    % ------------------------------------------------------------------
    %  Minimum is on LHS.  Done!  (Captures the dy=dz=0 case.)
    % ------------------------------------------------------------------
       sStat = STEP_LHS;
       cStep = 0;
    
    else
    % ------------------------------------------------------------------
    %  Minimizer must be to the right. The next breakpoint and its
    %  index is at the top of the heap; use heapdel to extract it and
    %  restore the the heap property to the remainder of the heap
    %  which will have length nBrk.
    % ------------------------------------------------------------------
       if using_octave                 % RHS of current interval
          [nBrk, rhs, iz, aBrk, iBrk] = heapdel(nBrk,aBrk,iBrk);
       else
          [nBrk, rhs, iz] = heapdel(nBrk,aBrk,iBrk);
       end
       if iz > n 
          iz = iz - n;                 % Recall that iz ranges over 1:2n
          tBreak = +1;                 % Breakpoint is an upper bound
       else
          tBreak = -1;                 % Breakpoint is a lower bound
       end
       lhs = tStep;                    % LHS of current interval
       maxStep = rhs - lhs;            % Max step to the next breakpoint
       
       if h < eps
           % This is a bit of a hack: a tiny h probably means that we're
           % already optimal (in BPdual).  This next line forces a quick
           % exit with a full steplength.
           sStat = STEP_MID;
%            % g is negative and q is nearly linear.
%            % Minimizer is either on RHS or at infinity.  Check.
%            if maxStep >= BigNum
%                sStat = STEP_INF;
%            else
%                sStat = STEP_RHS;
%                tStep = rhs;
%                cStep = maxStep;
%            end
       else
           % Minimum is on RHS or inside the interval. Check.
           cStep = -g / h;
           if cStep >= BigNum  &&  step < maxStep
               % Step is unbounded.  Exit with unbounded flag.
               sStat = STEP_INF;
           elseif cStep <= maxStep
               % Unconstrained minimizer found.  Exit.
               sStat = STEP_MID;
               tStep = tStep + cStep;
           else
               % Nondegenerate constrained minimizer on RHS.
               sStat = STEP_RHS;
               tStep = rhs;
               cStep = maxStep;
           end
       end
    end
        
    if debug
        if k==0
           fprintf(logH2,'k','iz','bound','type','cStep','tStep','bl',...
               'dz','z+tdz','bu','g','obj');
        end
        if   tBreak == +1, tBreakFlag = 'upp';
        else               tBreakFlag = 'low';
        end
        if   fixed(iz), typeFlag = 'fixd';
        else            typeFlag = 'free';
        end
        obj = funObj(tStep,tStep,dy);
        fprintf(logB2,k,iz,tBreakFlag,typeFlag,cStep,tStep,bl(iz),dz(iz),...
            z(iz)+tStep*(dz(iz)),bu(iz),g,obj)
    end
    
    % ------------------------------------------------------------------
    % Check exit conditions.
    % ------------------------------------------------------------------
    if sStat == STEP_RHS
        % We hit a breakpoint.  Keep going!
    else
        break
    end

    % ------------------------------------------------------------------
    % We'll now explore the next interval.  Update the gradient.
    % ------------------------------------------------------------------
    k = k + 1;
    diz  = dzt(iz);
    s.viz(k) = iz; % Save sequence of var indices.

    % Store current iterate info if this is the "best" breakpoint.
    if backtrack
        newangle = abs(angles(iz));
    else
        newangle = -inf;
    end
    if ~backtrack || (backtrack && newangle > s.aiz)
        s.k      = k;       % ...which breakpoint
        s.aiz    = newangle;
        s.tStep  = tStep;   % ...is the total steplength
        s.tBreak = tBreak;  % ...is it lower (-1) or upper (+1)
    end
    
    if diz == 0
        % Error: Can't have diz == 0 and a finite step.
        error('diz == 0 and STEP_RHS!!');
    elseif fixed(iz)
        g = g + rho * 2 * abs(diz);
        s.state(iz) = 2*sign(diz);
    elseif diz > 0
        g = g + rho * diz;
        if     s.state(iz) == -2
               s.state(iz) =   0;
        elseif s.state(iz) ==  0
               s.state(iz) =  +2;
        end
    elseif diz < 0
        g = g - rho * diz;
        if     s.state(iz) == +2
               s.state(iz) =   0;
        elseif s.state(iz) ==  0
               s.state(iz) =  -2;
        end
    else
        error('state(iz) = %i is inconsistent',iz,state(iz))
    end
    
    % Account for 2nd-order term.
    g = g + cStep * h;
    
end % while k <= nBrk
  
% ----------------------------------------------------------------------
% Three exit states are possible:
% 1. Unconstrained step and no breakpoints were encoutered.
%    No need to gather up the encountered breakpoints.
%    iz = 0
% 2. Unconstrained step and one or more breakpoints encountered.
%    Gather up the breakpoints stored in s.viz(1:k).
%    iz = index of last encountered breakpoint.
% 3. Constrained step.
%    Gather up breakpoints only up the 'best' breakpoint.
%    iz = index of 'best' encountered breakpoint.
% ----------------------------------------------------------------------
if k == 0
    tBreak = 0; iz = 0;
elseif sStat == STEP_MID
    tBreak = 0; iz = 0;
    state(s.viz(1:k)) = s.state(s.viz(1:k));   
else
    if backtrack
        k = s.k;
        iz = s.viz(k);
        tStep = s.tStep;
        tBreak = s.tBreak;
    end
    state(s.viz(1:k)) = s.state(s.viz(1:k));
    state(iz) = tBreak;
end
    
if debug, fprintf('\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRIVATE FUNCTIONS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function q = funObj(g,t,dy)
q = t*g + 0.5*t^2*norm(dy)^2;
end

