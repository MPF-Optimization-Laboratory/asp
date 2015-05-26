function [x,inform] = as_seqcs(MB, My, lambda, opts)
%AS_SEQCS  Sequential L1 for sparse signal recovery.
%
%   [X,INFORM] = AS_SEQCS(MB, MY, N) attempts to find the sparsest
%   representation of an N-by-1 signal y "sensed" via a linear measurement
%   operator.  If y is sparsely represented in terms of the n-by-N
%   basis/dictionary B, then
%
%     y = Bx  with  x  sparse.
%
%   The linear measurement operator M is used to measure y, and thus the
%   goal is to find a sparse solution to the system
%
%     M*B*x = M*y.
%
%   The input function MB returns products with the linear operator (M*B),
%   and the input function MY returns products M*y. These input functions
%   are described below.
%
%   AS_SEQCS(MB, My, N, OPTS) specifies options that can be set using
%   AS_SETPARMS.
%
%   The INFORM output argument is optional, and contains statistics on
%   the solution process.
%
%   Inputs
%   ------
%   MB is a function with the signature
%
%        A = MB(k)  which returns the compound linear operator A = M*B,
%   where
%        M is a  k-by-N linear operator (ie, the "measurement" matrix), and
%        B is an N-by-N linear operator (ie, the "sparsity basis" basis).
%
%   MY is a function with the signature
%
%        b = MY(k)  which returns the k-by-1 vector b = M*y.
%
%   N is the length of the unknown x.
%
%   OPTS is an options structure created using AS_SETPARMS.
%
%   See also AS_BPDN, AS_TOPY, AS_SETPARMS, BPDUAL.
%
%BPdual Toolbox
%Copyright 2010, Michael P. Friedlander and Michael A. Saunders
%http://www.cs.ubc.ca/labs/scl/bpdual

%$Id: as_seqcs.m 513 2010-04-20 16:43:12Z mpf $

REVISION = '$Revision: 310 $';
DATE     = '$Date: 2008-05-06 16:59:18 -0700 (Tue, 06 May 2008) $';
REVISION = REVISION(11:end-1);
DATE     = DATE(8:26);

% Check arguments
if nargin <  2 || isempty(MB) || isempty(My)
   error('At least 2 arguments needed');
end
if nargin < 3 || isempty(lambda), lambda = 0; end
if nargin < 4 || isempty(opts), opts = as_setparms; end

selftime = tic; % start the clock

% ---------------------------------------------------------------------
% Initialize local variables
% ---------------------------------------------------------------------
EXIT_CONVERGED       = 1;
EXIT_TOO_MANY_ITNS   = 2;
EXIT_BPERR           = 3;
EXIT_NO_RECOVERY     = 4;
exit_msg = {
   'Active set converged'
   'Too many iterations'
   'BPdual error'
   'No more observations available'
   };

eFlag    = false;
tol      = opts.tol_seqcs_convergence;      % stopping tolerance
k        = opts.num_additional_measurements;% additional measurements/itn
A        = MB(k);                           % Initially, k-by-n
m        = 0;                               % current no. of measurements
n        = size(A,2);
nprodA   = 0; nprodAt = 0;                  % no. of products with A
BPitns   = 0;                               % total no. of BPdual itns
BPtime   = 0;                               % total time spent in BPdual
BPopts   = as_setparms;                     % options for BPdual
BPopts.loglevel = opts.loglevel-1;
fid      = opts.fid;
[active,state,y,S,R] = deal([]);
itn      = 0;
itnMax   = floor(n/2);                      % rather arbitrary!
x        = zeros(n,1);

if lambda  < opts.lambdamin         % Don't allow a tiny value of lambda
   lambda  = opts.lambdamin;
   lamFlag = '!';                   % Flag that input lambda was ignored
else
   lamFlag = '';
end

%-----------------------------------------------------------------------
% Print log header.
%-----------------------------------------------------------------------
logB = ' %4i  %7i %7i %6.2f %11.4e %11.4e\n';
logH = ' %4s  %7s %7s %6s %11s %11s\n';

if opts.loglevel > 0
   fprintf(fid,'\n');
   fprintf(fid,' %s\n',repmat('=',1,80));
   fprintf(fid,' ASP: Sequential CS v.%s (%s)\n', REVISION, DATE);
   fprintf(fid,' %s\n',repmat('=',1,80));
   fprintf(fid,' %-20s: %8i %5s'    ,'No. coefficients'  ,n       ,'');
   fprintf(fid,' %-20s: %8.2e\n'    ,'Convergence tol'   ,tol        );
   fprintf(fid,' %-20s: %8i %5s'    ,'New measurements/itn',k,'');
   fprintf(fid,' %-20s: %8.2e%s\n'  ,'lambda'            ,lambda,lamFlag);
   fprintf(fid,'\n');
   fprintf(fid,logH,'Itn','BPitns','Active','(m/n)%','rNorm2','xNorm1');
end

%-----------------------------------------------------------------------
% Solve a sequence of BP problems
%-----------------------------------------------------------------------
while true
   
   itn = itn + 1;
   
   mpk = min(m+k,n);
   if itn > 1      % Already got operator above (needed to infer n).
      A = MB(mpk); % New operator with an add'l k rows: now (m+k)-by-n
   end
   b = My(mpk);    % New RHS vector: now (m+k)-by-1
   
   % For each new row in A, update Q-less QR factorization of S.
   if ~isempty(active)
      y(m+1:mpk,:) = 0; % this indexing forces y to expand columnwise
      z = zeros(mpk,1);
      for i = m+1:mpk
         z(i)   = 1;
         a      = (A'*z)';  nprodAt = nprodAt + 1;
         z(i)   = 0;
         s      = a(active);
         S(i,:) = s;
         R      = QRaddrow(R,s);
      end
   end
   m = mpk;
   
   % Solve the current subproblem
   [active,state,xx,y,S,R,inform] = ...
      BPdual(A,b,-1,1,lambda,active,state,y,S,R,BPopts);
   
   % Grab BPdual stats and parameters
   xNorm   = inform.xNorm;
   rNorm   = inform.rNorm;
   nact    = length(active);
   BPtime  = BPtime + inform.time;
   BPitns  = BPitns + inform.itns;
   nprodA  = nprodA  + inform.nprodA;
   nprodAt = nprodAt + inform.nprodAt;
   
   % Save the old solution
   xold = x;
   
   % BPdual's solution x is short. Make it full length.
   x = zeros(n,1);
   x(active) = xx;
   
   % Check if BPdual failed or max number of iterations.
   if inform.itns == 0
      % Presume that additional observations aren't changing the solution.
      eFlag = EXIT_CONVERGED;
   elseif norm(x-xold,inf) < max(1,norm(xold,inf))*tol
      eFlag = EXIT_CONVERGED;
   elseif inform.stat == 0 || inform.stat == 6
      % No BPdual errors. Keep going
      eFlag = 0;
   elseif inform.stat
      % Unhandles BPdual error.
      eFlag = EXIT_BPERR;
   elseif itn >= itnMax
      % Ran out of iterations.
      eFlag = EXIT_TOO_MANY_ITNS;
   elseif m == n
      % No more observations available! Not converged.
      eFlag = EXIT_NO_RECOVERY;
   end
   
   % Print to log
   if opts.loglevel > 0
      fprintf(fid,logB,itn,BPitns,nact,m/n*100,rNorm,xNorm);
   end
   
   % Act on exit conditions
   if eFlag, break, end
   
end % while


inform.lam     = lambda;
inform.y       = y;
inform.active  = active;
inform.state   = state;
inform.S       = S;
inform.R       = R;
inform.itns    = itn;
inform.BPitns  = BPitns;
inform.nprodA  = nprodA;
inform.nprodAt = nprodAt;
inform.BPtime  = BPtime;
inform.nobs    = m;
inform.time    = toc(selftime);

if opts.loglevel > 0
   fprintf(fid,'\n EXIT %s -- %s\n',mfilename,exit_msg{eFlag});
   fprintf(fid,' %-20s: %8i %5s','No. significant nnz',sparsity(xx),'');
   fprintf(fid,' %-20s: %8i\n','Products with A',nprodA);
   fprintf(fid,' %-20s: %8.1e %5s','Solution time (sec)',inform.time,'');
   fprintf(fid,' %-20s: %8i\n','Products with At',nprodAt);
   fprintf(fid,' %-20s: %8.1e %5s','BPdual time (sec)',BPtime,'');
   fprintf(fid,' %-20s: %8i\n','BPdual iterations',BPitns);
   fprintf(fid,' %-20s  %8.1e %5s','Final residual',rNorm,'');
   fprintf(fid,'\n');
end

%---------------------------------------------------------------------
end % function as_rwbp.
%---------------------------------------------------------------------
