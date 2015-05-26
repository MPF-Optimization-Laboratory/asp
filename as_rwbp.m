function [x,inform] = as_rwbp(A, b, opts)
%AS_RWBP  Reweighted L1 for sparse Ax=b.
%
%   [X,INFORM] = AS_RWPB(A,B) attempts to find the sparsest solution
%   to the linear system
%
%      Ax = b.
%
%   AS_RWPB(A,B,OPTS) specifies options that can be set using
%   AS_SETPARMS.
%
%   The INFORM output argument is optional, and contains statistics on
%   the solution process.
%
%   Inputs
%   A       is an m-by-n matrix, explicit or an operator.
%           If A is a function, then it must have the signature
%
%           y = A(x,mode)   if mode == 1 then y = A x  (y is m-by-1);
%                           if mode == 2 then y = A'x  (y is n-by-1).
%
%   B       is an m-vector.
%   OPTS    is an options structure created using AS_SETPARMS.
%
%   Example
%   m = 600; n = 2560; k = 20;    % No. of rows, columns, and nonzeros
%   p = randperm(n); p = p(1:k);  % Position of nonzeros in x
%   x = zeros(n,1);               % Generate sparse solution
%   x(p) = randn(k,1);
%   A = randn(m,n);               % Gaussian m-by-n ensemble
%   b = A*x;                      % Compute the RHS vector
%   [x,inform] = as_rwpb(A,b);    % Solve the basis pursuit problem
%
%   See also AS_BPDN, AS_TOPY, AS_SETPARMS, BPDUAL.
%
%BPdual Toolbox
%Copyright 2008, Michael P. Friedlander and Michael A. Saunders
%http://www.cs.ubc.ca/labs/scl/bpdual

%$Id: as_rwbp.m 451 2009-03-23 15:13:25Z mpf $

REVISION = '$Revision: 310 $';
DATE     = '$Date: 2008-05-06 16:59:18 -0700 (Tue, 06 May 2008) $';
REVISION = REVISION(11:end-1);
DATE     = DATE(8:26);

% Check arguments
if nargin <  2, error('At least 2 arguments needed'); end
if nargin <  3, opts = []; end

tic;             % start the clock
[m,n] = size(A); % Grab the number of rows and columns.

% Default options
if isempty(opts)
    opts = as_setparms;
end

% BPdual gets its own set of options.
BPopts.loglevel = opts.loglevel-1;
BPopts.elastic = false;       % initially fire up BPdual in inelastic mode
BPopts = as_setparms(BPopts); % Set defaults for remainder of options

% Initialize local variables
itn      = 0;
itnMax   = opts.rwitns;
rwfid    = opts.rwfid;
rwtol    = opts.rwtol;
feaTol   = opts.rwfeaTol;
w        = ones(n,1);      % Weights: bl = -w; bu = w;
[active,state,y,S,R] = deal([]);
if opts.rwtopy
    z   = A'*b;
    lam = norm(z,inf) / 2;
else
    lam = 0;
end

% Exit flags
EXIT_NOERR  = 0;
EXIT_BPERR  = 1;
EXIT_ITNS   = 2;
EXIT_SPARK  = 3;
eFlag = EXIT_NOERR;

% LSQR parameters
atol = feaTol; btol = feaTol; conlim = 1e+6; damp = 0;
CGmaxits = m; CGshow = 0;

% Statistics
BPitns  = 0; BPtime  = 0; CGitns = 0;
nprodA  = 0; nprodAt = 1;  % Already one prod with A' above

%-----------------------------------------------------------------------
% Print log header.
%-----------------------------------------------------------------------
logB = ' %4i  %7i %7i %11.4e %11.4e %11.4e %11.4e %11.4e\n';
logH = ' %4s  %7s %7s %11s %11s %11s %11s %11s\n';

if opts.loglevel > 0
   fprintf(rwfid,'\n');
   fprintf(rwfid,' %s\n',repmat('=',1,80));
   fprintf(rwfid,' ASP: Rewieghted BP v.%s (%s)\n', REVISION, DATE);
   fprintf(rwfid,' %s\n',repmat('=',1,80));
   fprintf(rwfid,' %-20s: %8i %5s'    ,'No. rows'          ,m       ,'');
   fprintf(rwfid,' %-20s: %8.2e\n'    ,'Reweight tol'      ,rwtol      );
   fprintf(rwfid,' %-20s: %8i %5s'    ,'No. columns'       ,n       ,'');
   fprintf(rwfid,' %-20s: %8i\n'      ,'Maximum iterations',itnMax     );
   fprintf(rwfid,'\n');
   fprintf(rwfid,logH,'Itn','BPitns','Active','Weight min','Weight max',...
       'lambda','rNorm2','xNorm1');
end

%-----------------------------------------------------------------------
% Solve a sequence of weighted BP problems
%-----------------------------------------------------------------------
while true

    itn = itn + 1;
    
    [active,state,xx,y,S,R,inform] = ...
        BPdual(A,b,-w,w,lam,active,state,y,S,R,BPopts);

    % Grab BPdual stats and parameters
    xNorm   = inform.xNorm;
    rNorm   = inform.rNorm;
    nact    = length(active);
    wmin    = min(w); wmax = max(w);
    rho     = inform.rho;
    BPitns  = BPitns  + inform.itns;
    BPtime  = BPtime  + inform.time;
    nprodA  = nprodA  + inform.nprodA;
    nprodAt = nprodAt + inform.nprodAt;
    
    % Print to log
    if opts.loglevel > 0
       if opts.loglevel > 1, fprintf(rwfid,'\n'); end
       fprintf(rwfid,logB,itn,BPitns,nact,wmin,wmax,lam,rNorm,xNorm);
    end
    
    % Check if BPdual failed or max number of iterations.
    if inform.stat ~= 0 && inform.stat ~= 6
       eFlag = EXIT_BPERR;
    elseif itn >= itnMax
       eFlag = EXIT_ITNS;
    end
    
    % If A has full rank then the current support is the sparsest
    % possible if it has fewer than m/2 elements, and the following
    % least-squares problem has zero residual.
%   if 0
    if nact < m/2 && (rNorm < sqrt(feaTol) || lam < sqrt(feaTol))
%   if rNorm*lam < 1e-6 && nact < m/2
%   if nact < m/2
       Aind = @(x,mode)aprodind(mode,x,A,m,n,active);
       [xxCG,istop,CGitnk,CGrnorm] = ...
           lsqr(m,nact,Aind,b,damp,atol,btol,conlim,CGmaxits,CGshow);
       CGitns     = CGitns  + CGitnk;
       nprodA     = nprodA  + CGitnk;
       nprodAt    = nprodAt + CGitnk + 1;
       if istop == 1 % The system A(:,active)xx=b is compatible.
           xx    = xxCG;
           rNorm = CGrnorm;
           eFlag = EXIT_SPARK;
       end
    end

    % Scatter solution into full vector. (Needed to set weights below.)
    x = zeros(n,1);
    x(active) = xx;

    if eFlag, break, end              % Act on error flags.
    
    if opts.rwtopy, lam = lam/2; end  % Smaller lambda    
    w = 1./(abs(x) + rwtol);          % Set new weights
    BPopts.rho = max(xNorm,rho);   % Penalty param for next elastic phase
end

switch eFlag
    case EXIT_NOERR
        exit_msg = sprintf('No BPdual lerrors\n');
    case EXIT_SPARK
        exit_msg = sprintf('Sparest solution found\n');
    case EXIT_ITNS
        exit_msg = sprintf('Too many iterations');
    case EXIT_BPERR
        exit_msg = sprintf('BPdual error exit %i -- %s\n',...
            inform.stat, inform.exitmsg);
end

inform.lam     = lam;
inform.y       = y;
inform.active  = active;
inform.state   = state;
inform.S       = S;
inform.R       = R;
inform.rwitns  = itn;
inform.BPitns  = BPitns;
inform.nprodA  = nprodA;
inform.nprodAt = nprodAt;
inform.rwtime  = toc;
inform.BPtime  = BPtime;
inform.error   = exit_msg;

if opts.loglevel > 0
   fprintf(rwfid,'\n EXIT as_rwbp -- %s\n',exit_msg);
   fprintf(rwfid,' %-20s: %8i %5s','No. significant nnz',sparsity(xx),'');
   fprintf(rwfid,' %-20s: %8i\n','Products with A',nprodA);
   fprintf(rwfid,' %-20s: %8.1e %5s','Solution time (sec)',inform.rwtime,'');
   fprintf(rwfid,' %-20s: %8i\n','Products with At',nprodAt);
   fprintf(rwfid,' %-20s: %8.1e %5s','BPdual time (sec)',BPtime,'');
   fprintf(rwfid,' %-20s: %8i\n','BPdual iterations',BPitns);
   fprintf(rwfid,' %-20s  %8.1e %5s','Final residual',rNorm,'');
   fprintf(rwfid,' %-20s: %8i\n','LSQR iterations',CGitns);
   fprintf(rwfid,'\n');
end

%---------------------------------------------------------------------
end % function as_rwbp.
%---------------------------------------------------------------------
