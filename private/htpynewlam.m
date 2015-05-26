function [x,dy,dz,step,lambda,p] = ...
        htpynewlam(active,state,A,R,S,x,y,s1,s2,lambda,lamFinal)
%htpynewlam  finds a new lambda for the homotopy method.

% 06 Jul 09: Allow for general bounds.
% 04 Aug 09: Xiangrui Meng found that we can get alfamin < 0
%            when BPdual has lambdamin = 0.  (Possibly other
%            times too -- we aren't sure if s1,s2 >= 0 strictly.
%            Not yet fixed.
% 20 Aug 09: If "next" lambda lies below lamFinal, compute step so that
%            lambda=lamFinal.

% $Id: htpynewlam.m 490 2009-08-20 17:20:47Z mpf $

p    = [];          % will be index of new active constraint
alfa =   inf(4,1);
k    = zeros(4,1);

% Compute the search directions; see (6).
[dx,dy] = csne(R,S,y);
dz      = A'*dy;

% Find the largest allowable step along dx without sign changes.
I1 = find( state(active)==+1  &  dx < -eps );    % x  pos and decreasing
I2 = find( state(active)==-1  &  dx >  eps );    % x  neg and increasing
if any(I1)
   [alfa(1),i] = min(-x(I1)./dx(I1)); i1=I1(i); k(1)=active(i1);
end
if any(I2)
   [alfa(2),i] = min(-x(I2)./dx(I2)); i2=I2(i); k(2)=active(i2);
end

% Find the largest allowable step along dy without violating a constraint.
free = abs(state) ~= 1;
r1 = dz + s1;
r2 = dz - s2;
I3 = find( dz < -eps  &  free );
I4 = find( dz >  eps  &  free );
if any(I3)
   [alfa(3),i] = min( lambda*s1(I3)./r1(I3)); k(3) = I3(i);
end
if any(I4)
   [alfa(4),i] = min(-lambda*s2(I4)./r2(I4)); k(4) = I4(i);
end
  
% Reduce lambda and update variables
[alfamin, kk] = min(alfa);

% Lambda should never increase, ie, alphamin should always be nonnegative.
if alfamin < 0
   if alfamin > -1e-15
      % Don't panic -- it seems to be a degenerate step and roundoff.
      % Explicitly set the step to 0.
      alfamin = 0;
   else
      error('%s: change in lambda is negative: %s -- call for help!',...
         mfilename,alfamin)
   end
end   

alfamax = lambda - lamFinal;
dlambda = min( alfamin, alfamax );
lamNew  = lambda - dlambda;

% Decide if we need to force an active-set change
if lamNew > lamFinal   % Still haven't reached the limit on lambda
   force  = true;      % We know an active-set change will occur. Force it
else
   force  = false;     % No active-set change will occur. Don't push it
end

% Update the variables.
x = x + dlambda*dx;
if norm(dy) <= (1+norm(y))*sqrt(eps)
   step = 0;
else
   step = dlambda/lamNew;
end

% If a variable is blocking, then unambiguously force the sign.
% If a constraint is blocking, then augment the step to force a hit bound.
if force
   if     kk == 1       % x will go neg
      x(i1) = -1;
   elseif kk == 2       % x will go pos
      x(i2) = +1;
   elseif kk == 3       % A new lower-bound becomes active
      p = -k(kk);      % ...this is its index.
   else
      p =  k(kk);
   end
end

% Set lambda to its new value.
lambda = lamNew;
