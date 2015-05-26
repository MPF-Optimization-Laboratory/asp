function [active,state,x,r,S,R,inform] = as_omp(A,b,lambda,active,state,r,S,R,opts)
%AS_OMP  Orthogonal matching pursuit for sparse Ax=b.
%
% OMP applies the orthogonal matching pursuit (OMP) algorithmn to
% estimate a sparse solution of the underdetermined system Ax=b.
%
% Solve
% [active,state,x,r,S,R,inform] = as_omp(A,b,lambda);
%
% Restart (with some other lambda):
% [active,state,x,r,S,R,inform] = as_omp(A,b,lambda,active,state,r,S,R);
%
% INPUTS
% ======
% A          is a function handle representing an m x n matrix "A"
%            such that
%              A(x,1) returns A*x
%              A(x,2) returns A'*x
% b          is an m-vector.
% lambda     is the smallest allowed value of ||A'r||_inf.
% active,state,r,S,R  may be empty
%            or output from OMP with some previous value of lambda.
%            
% OUTPUT
% active     is an nact-vector of indices j with nact = length(active),
%            listing which of the "constraints"  -lam e <= A'r <= lam e
%            is active.
% state      is an n-vector giving the state of each dual constraint.
%            state(j) = -lambda  if         Aj'*y =  lam
%                     = +lambda  if         Aj'*y = -lam
%                     =  0       if   lam < Aj'*y <  lam
% x          is an nact-vector of solution values.
%            The n-vector xx = zeros(n,1), xx(active) = x.
% r          is an m-vector with the current residual.
% S          is the submatrix A(:,active).
% R          is a Cholesky factorization of S.
% inform     is a structure with the following components:
%            .time    total solution time (seconds)
%            .stat    = 0  solution is optimal
%                     = 1  b = 0                ---> x,y = 0
%                     = 2  lambda > ||A'b||_inf ---> x,y = 0
%                     = 3  dual infeasible point
%                     = 4  too many iterations
%                     = 5  current S is rank-deficient.
%                     = 6  r has gone infeasible.

% 09 Feb 2008: First version of OMP.
%              Michael Friedlander and Michael Saunders.
%              Least-squares subproblems solved by x = S\b.
%
%----------------------------------------------------------------------
  REVISION = '$Revision: 223 $';
  DATE     = '$Date: 2008-02-08 19:33:59 -0800 (Fri, 08 Feb 2008) $';
  REVISION = REVISION(11:end-1);
  DATE     = DATE(8:26);

%----------------------------------------------------------------------
% Check input arguments.
%----------------------------------------------------------------------
  if nargin < 2 || nargin > 9
      error('At least two arguments are required');
  end
  if nargin < 8
      opts = [];
  end
  if nargin < 4
      active = [];  % Forces a cold start.
  end
  if nargin < 3 || isempty(lambda) || lambda <= 0
      lambda = sqrt(eps);
  end

%-----------------------------------------------------------------------
% Start the clock and size up the problem.
%-----------------------------------------------------------------------
  tic;
  z        = A'*b;
  m        = length(b);
  n        = length(z);
  nprodA   = 0;
  nprodAt  = 1;

%----------------------------------------------------------------------
% Grab input options and set defaults where needed. 
%----------------------------------------------------------------------
  if isempty(opts)
      opts = as_setparms;
  end
  itnMax   = opts.ompitns*m;
  opttol   = opts.opttol;
  
  printf(opts,'\n');
  printf(opts,' %s\n',repmat('=',1,80));
  printf(opts,' OMP  v.%s (%s)\n', REVISION, DATE);
  printf(opts,' %s\n',repmat('=',1,80));
  printf(opts,' %-20s: %8i %5s'    ,'No. rows'          ,m       ,'');
  printf(opts,' %-20s: %8.2e\n'    ,'lambda'            ,lambda     );
  printf(opts,' %-20s: %8i %5s'    ,'No. columns'       ,n       ,'');
  printf(opts,' %-20s: %8.2e\n'    ,'Optimality tol'    ,opttol     );
  printf(opts,' %-20s: %8i %5s'    ,'Maximum iterations',itnMax  ,'');
  printf(opts,'\n');

%----------------------------------------------------------------------
% Initialize local variables.
%----------------------------------------------------------------------
  EXIT_UNDEFINED       = 0;    % Exit flags.
  EXIT_OPTIMAL         = 1;
  EXIT_TOO_MANY_ITNS   = 2;
  EXIT_SINGULAR_LS     = 3;
  EXIT_LAMBDA          = 4;
  EXIT_RHS_ZERO        = 5;
  EXIT_UNCONSTRAINED   = 6;
  exit_msg = {
      'Optimal solution found'
      'Too many iterations'
      'Singular least-squares subproblem'
      'Reached minimum value of lambda'
      'b = 0. The solution is x = 0'
      'Unconstrained solution  r = b  is optimal'
             };

  itn     = 0;
  eFlag   = EXIT_UNDEFINED;
  x       = zeros(0,1);
  zerovec = zeros(n,1);
  p       = 0;

% Quick exit if the RHS is zero.
  if norm(b,inf) == 0
     r      = zeros(m,1);
     eFlag  = EXIT_RHS_ZERO;
  end

% Solution is unconstrained for lambda large.
  zmax = norm(z,inf);
  if ~eFlag && zmax < lambda
     r      = b;
     eFlag  = EXIT_UNCONSTRAINED;
  end

% Reset the active set indices if this is a cold start or there's an error.
  if eFlag || isempty(active)
      active = zeros(0,1);
      state  = zeros(n,1);
      S      = zeros(m,0);
      R      = zeros(0,0);
  end

  logB = ' %4i  %8i %12.5e %12.5e %12.5e\n';
  logH = ' %4s  %8s %12s %12s %12s\n';
  printf(opts,'\n');
  printf(opts,logH,'Itn','Var','lambda','rNorm','xNorm')

%-----------------------------------------------------------------------
% Main loop.
%-----------------------------------------------------------------------
  while 1
    %-------------------------------------------------------------------
    % Compute dual obj gradient g, search direction dy, and residual r.
    %-------------------------------------------------------------------
    if itn == 0
       x    = [];
       r    = b;
       z    = A'*r;  nprodAt = nprodAt+1;
       zmax = norm(z,inf);
    else
       x = csne(R,S,b);
       if norm(x,inf) > 1e+12
          eFlag = EXIT_SINGULAR_LS;
          break
       end
       Sx   = S*x;
       r    = b - Sx;
    end
    
    rNorm = norm(r,2);
    xNorm = norm(x,1);
    
    printf(opts,logB,itn,p,zmax,rNorm,xNorm);

    %---------------------------------
    % Check exit conditions.
    %---------------------------------
    if eFlag
        % Already set. Don't test the other exits.
    elseif zmax <= lambda
        eFlag = EXIT_LAMBDA;
    elseif rNorm <= opttol
        eFlag = EXIT_OPTIMAL;
    elseif itn >= itnMax
        eFlag = EXIT_TOO_MANY_ITNS;
    end
    if eFlag, break, end
    
    %-------------------------------------------
    % New iteration starts here.
    %-------------------------------------------
    itn    = itn + 1;

    %---------------------------------------------
    % Find step to the nearest inactive constraint
    %---------------------------------------------
    z        = A'*r;         nprodAt = nprodAt+1;
    [zmax,p] = max(abs(z));
    if z < 0
        state(p) = -1;
    else
        state(p) = +1;
    end
    zerovec(p) = 1;           % Extract a = A(:,p)
    a          = A*zerovec;   nprodA = nprodA+1;
    zerovec(p) = 0;
    R          = QRaddcol(S,R,a);% Update R
    S          = [ S       a ];  % Expand S, active
    active     = [ active; p ];

  end % while 1
%-----------------------------------------------------
% end main loop
%-----------------------------------------------------

  inform.itns    = itn;                   % Gather exit data.
  inform.time    = toc;
  inform.stat    = eFlag - 1;
  inform.exitmsg = exit_msg{eFlag};
  inform.nprodA  = nprodA;
  inform.nprodAt = nprodAt;
  
  printf(opts,'\n EXIT OMP -- %s\n\n',exit_msg{eFlag});
  printf(opts,' %-20s: %8i %5s','No. significant nnz',sparsity(x),'');
  printf(opts,' %-20s: %8i\n','Products with A',nprodA);
  printf(opts,' %-20s: %8.1e %5s','Solution time (sec)',inform.time,'');
  printf(opts,' %-20s: %8i\n','Products with At',nprodAt);

%---------------------------------------------------------------------
end % function omp.
%---------------------------------------------------------------------

function printf(opts,varargin)
  if opts.loglevel > 0
     fprintf(opts.fid,varargin{:});
  end
end % function printf
