function [H,noup] = lbfgsupdate(H, step, p, g1, g2)
%lbfgsupdate  updates the L-BFGS matrix.
%
%   [H,noup] = lbfgsupdate(H,step,p,g1,g2) updates H.  The flag noup
%   indicates if the update was skipped.
%
%   See also lbfgsadd, lbfgsdel, lbfgshprod, lbfgsinit.

% 02 Sept 2008: Continue to "age" the vectors even if an update is skipped.
%               This means that the pointers will get advanced and that the
%               oldest vectors won't participate in the next Hprod.

% Storage scheme: jOld points to the location with the oldest
% vectors. The next location has the newest vectors.  For examples, if
% jMax = 5, then
%
% |---|---|---|---|---|
%       ^   ^
%       |   |
%      old new
%   4   5   1   2   3    1=newest, 2=2nd newest,...,5=oldest
%
% $Id: lbfgsupdate.m 489 2009-08-20 16:43:09Z mpf $

% Get size of storage and pointer to the last updated vectors.
jMax = H.jMax;
jOld = H.jOld;

gtp1 = g1'*p;
gtp2 = g2'*p;
noup = gtp2 <= 0.91*gtp1;

if noup
   H.r(jOld) = 0;    % Causes the corresponding vectors to not participate
                     % when computing subsequent products with H. 
else
   y = g2 - g1;

   % Replace oldest vectors in S and Y.
   H.S(:,jOld) = step*p;
   H.Y(:,jOld) = y;
   H.r(  jOld) = 1 / (step*(gtp2 - gtp1));

   % Update diagonal
   gamma   = 1/(H.r(jOld)*norm(y)^2); % = s'y/y'y
   H.H0(:) = gamma;
end

% Advance pointers: the newest vectors are now where the oldest used to be.
% The oldest vectors are always one spot "behind" the newest.
H.jNew = jOld;
if jOld == 1
   H.jOld = jMax;
else
   H.jOld = jOld - 1;
end

% ---------------------
% DEBUG - DEBUG - DEBUG
% ---------------------

%disp(H.rank)
%g0 = randn(size(H.S,1),1);
%g1 = lbfgsbprod(H,g0);
%g2 = lbfgshprod(H,g1);
%disp(norm(g2-g0))
%plot(g0); hold on; plot(g2,'.'); hold off; pause;
