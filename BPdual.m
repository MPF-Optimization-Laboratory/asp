function [active,state,x,y,S,R,inform] = BPdual(...
   A,b,bl,bu,lambdain,active,state,y,S,R,opts)
% BPdual solves the optimization problem
%
%    DP:       minimize_y   - b'y  +  1/2 lambda y'y
%              subject to   bl <= A'y <= bu
%
% using given A,b,bl,bu,lambda.  When bl = -e, bu = e = ones(n,1),
% DP is the dual Basis Pursuit problem
%
%    BPdual:   maximize_y   b'y  -  1/2 lambda y'y
%              subject to   ||A'y||_inf  <=  1.
%
% General call to BPdual:
% [active,state,x,y,S,R,info] = BPdual(A,b,bl,bu,lambda,active,state,y,S,R,opts)
%
% EXAMPLE
% =======
% % Generate data:
% [A,b,lambda] = gendata(1);
% % Solve BPdual:
% [active,state,x,y,S,info] = BPdual(A,b,-1,1,lambda);
% % Resolve with a different lambda
% lambda = lambda / 10;
% [active,state,x,y,S,R,info] = BPdual(A,b,-1,1,lambda,active,state,y,S,R);
%
% INPUT
% =====
% A          is an m-by-n explicit matrix or an operator.
% b          is an m-vector.
% lambda     is a nonnegative scalar.
% bl, bu     are n-vectors. Optionally, they can be scalar valued,
%            and in that case they are assumed to be constant n-vectors.
% active,state,y,S,R  may be empty (i.e., equal to [])
%            or output from BPdual with some previous value of lambda.
%            
% OUTPUT
% ======
% active     is an nact-vector of indices j with nact = length(active),
%            listing which of the constraints bl(j) <= A(:,j)'y <= bu(j)
%            is active.
% state      is an n-vector giving the state of each dual constraint.
%            state(j) = -2  if         Aj'*y < bl(j)  infeas below
%                     = -1  if         Aj'*y = bl(j)  active below
%                     = +1  if         Aj'*y = bu(j)  active above
%                     = +2  if         Aj'*y > bu(j)  infeas above
%                     =  0  if bl(j) < Aj'*y < bu(j)  inactive
% x          is an nact-vector of primal solution values.
%            The n-vector xx = zeros(n,1), xx(active) = x
%            solves the primal BP problem (see below).
% y          is an m-vector of BPdual solution values.
% S          is the submatrix A(:,active).
% R          is the Cholesky factorization of S.
% info       is a structure with the following components:
%            .time    total solution time (seconds)
%            .stat    = 0  solution is optimal
%                     = 1  too many iterations
%                     = 2  current S is rank-deficient
%                     = 3  dual infeasible point
%            .itn     number of iterations
%            .exitmsg exit message
%            .nprodA  number of products with A
%            .nprodAt number of products with A'
%            .rho     final elastic penalty parameter (may not apply)
%            .rNorm   2-norm of the final residual r = b - A(:,active)*x
%            .xNorm   1-norm of the primal solution
%            .opts    options used
%            .numtrim number of constraints trimmed from final working set

%BPdual Toolbox
%Copyright 2008, Michael P. Friedlander and Michael A. Saunders
%http://www.cs.ubc.ca/labs/scl/bpdual

% The primal BP problem may be written
%
%    BP:        min_{x,y} ||x||_1 +  1/2 lambda||y||_2^2
%               st         Ax     +      lambda  y  =  b.
%
% If lambda > 0, we have the residual vector r = lambda y = b - Ax,
% and BP is equivalent to the original BP-denoising problem
%
%    BPDN:      min_{x,y} lambda ||x||_1 +  1/2 ||r||_2^2
%               st         Ax     +      r  =  b.
%
%
% 09 Jun 2007: First version of BPdual.
%              Michael Friedlander and Michael Saunders.
%              A can be dense or sparse (but not an operator).
%              S = A(:,active) holds active columns.
%              Least-squares subproblems solved by x = S\g.
% 29 Aug 2007: "A" is now an operator.
% 09 Nov 2007: Generalize the l1 dual constraints from -e <= A'y <= e
%              to bl <= A'y <= bu.
%
% $Id: BPdual.m 514 2010-04-20 16:47:44Z mpf $

  REVISION = '$Revision: 514 $';
  DATE     = '$Date: 2010-04-20 09:47:44 -0700 (Tue, 20 Apr 2010) $';
  REVISION = REVISION(11:end-1);
  DATE     = DATE(8:26);
  selftime = tic;

%----------------------------------------------------------------------
% Check input arguments.
%----------------------------------------------------------------------
  if nargin < 5 || nargin > 11
     error('Wrong number of arguments.');
  end
  if nargin < 6 ||isempty(active)||isempty(y)||isempty(S)||isempty(R)
     coldstart = true;
  else
     coldstart = false;
  end
  if nargin < 11 || isempty(opts)
     opts = as_setparms;
  end

  % 11 Jun 09: Right now BPdual only allows coldstarts (ie, y=0) with the
  % homotopy method. It should allow for homotopy from some previous call
  % to BPdual. FIX!
  if coldstart || opts.homotopy
     z = A'*b;      % Normally, z = A'y. This is the one exception.
  else
     z = A'*y;
  end
  m        = length(b);
  n        = length(z);
  nprodA   = 0;
  nprodAt  = 1;

  % Accept scalars for the upper and lower bounds. But it simplifies
  % code if we immediately turn them into vectors.
  if length(bl) == 1, bl = repmat(bl,n,1); end
  if length(bu) == 1, bu = repmat(bu,n,1); end

%----------------------------------------------------------------------
% Grab input options and set defaults where needed. 
%----------------------------------------------------------------------
  fid       = opts.fid;
  itnMax    = opts.majiters*max(m,n);
  elastic   = opts.elastic;
  homotopy  = opts.homotopy;
  featol    = opts.featol;
  gaptol    = opts.gaptol;
  opttol    = opts.opttol;
  pivtol    = opts.pivtol;
  tietol    = featol + 1e-8;           % Perturbation to break ties
  rhomax    = opts.rhomax;
  rho       = opts.rho;
  actmax    = opts.actmax;
  loglevel  = opts.loglevel;
  lambdamin = opts.lambdamin;
  trim      = opts.trim;
  callback  = isa(opts.callback,'function_handle');
  if lambdain < lambdamin              % Don't allow a tiny value of lambda
      lambda  = lambdamin;
      lamFlag = '!';                   % Flag that input lambda was ignored
  else
      lambda  =  lambdain;
      lamFlag = '';
  end
  if homotopy
      elastic  = false;                % Don't allow elatic mode
      gaptol   = 0;                    % Ignore gaptol: each itn is optimal
      lamFinal = lambda;               % Lambda might have been adjusted
      lambda   = norm(z,inf);          % Currently, z = A'b
  end

%-----------------------------------------------------------------------
% Print log header.
%-----------------------------------------------------------------------
  logB = ' %4i  %8.1e %8s %7i %11.4e %11.4e %11.4e %7.1e  %7.1e %4i %7.1e';
  logH = ' %4s  %8s %9s %7s %11s %11s %11s %7s  %7s %4s %7s';
  if loglevel > 0
      fprintf(fid,'\n');
      fprintf(fid,' %s\n',repmat('=',1,80));
      fprintf(fid,' BPdual  v.%s (%s)\n', REVISION, DATE);
      fprintf(fid,' %s\n',repmat('=',1,80));
      fprintf(fid,' %-20s: %8i %5s'    ,'No. rows'          ,m       ,'');
      fprintf(fid,' %-20s: %8.2e%s\n'  ,'lambda'            ,lambda,lamFlag);
      fprintf(fid,' %-20s: %8i %5s'    ,'No. columns'       ,n       ,'');
      fprintf(fid,' %-20s: %8.2e\n'    ,'Optimality tol'    ,opttol     );
      fprintf(fid,' %-20s: %8i %5s'    ,'Maximum iterations',itnMax  ,'');
      fprintf(fid,' %-20s: %8.2e\n'    ,'Duality tol'       ,gaptol     );
      fprintf(fid,' %-20s: %8i %5s'    ,'Support trimming'  ,trim,    '');
      fprintf(fid,' %-20s: %8.2e\n'    ,'Pivot tol'         ,pivtol     );
      if elastic
      fprintf(fid,' %-20s: %8.2e %5s'  ,'Initial penalty'   ,rho     ,'');
      fprintf(fid,' %-20s: %8.2e\n'    ,'Maximum penalty'   ,rhomax     );
      end
  end

%----------------------------------------------------------------------
% Initialize local variables.
%----------------------------------------------------------------------
  EXIT_OPTIMAL         = 1;
  EXIT_TOO_MANY_ITNS   = 2;
  EXIT_SINGULAR_LS     = 3;
  EXIT_INFEASIBLE      = 4;
  EXIT_REQUESTED       = 5;
  EXIT_ACTMAX          = 6;
  EXIT_SMALL_DGAP      = 7;
  exit_msg = {
      'Optimal solution found -- full Newton step'
      'Too many iterations'
      'Singular least-squares subproblem'
      'Reached dual infeasible point'
      'User requested exit'
      'Max no. of active constraints reached'
      'Optimal solution found -- small duality gap'
             };

  itn       = 0;
  step      = 0;
  p         = [];                        % index of added constraint
  q         = [];                        % index of deleted constraint
  svar      = '';                        % ...and their string value.
  eFlag     = false;
  zerovec   = zeros(n,1);
  nBrks     = 0;
  newInfeas = true;
  phase1    = false;
  numtrim   = 0;
% cNorms = colnorms(A,m,n);              % Currently not needed

%-----------------------------------------------------------------------
% Cold/warm-start initialization.
%-----------------------------------------------------------------------
  if coldstart
     active = zeros(0,1);
     state  = zeros(n,1);
     if issparse(A), S = sparse(m,0);
     else            S =  zeros(m,0);
     end
     R = zeros(0,0);
     x = zeros(0,1);
     if homotopy
         y  = b / lambda;
         z  = z / lambda;              % This sets z = A'y
     else
         y  = zeros(m,1);
         z  = zeros(n,1);
     end
  else                                 % This is a warmstart
     % Trim active-set components that correspond to huge bounds, and then
     % find the min-norm change to y that puts the active constraints
     % exactly on their bounds.
     g    = b - lambda*y;              % Compute steepest-descent dir
     [active,state,S,R] = triminf(active,state,S,R,bl,bu,g,opts);
     nact = length(active);
     x    = zeros(nact,1);
     y    = restorefeas(y,active,state,S,R,bl,bu);     
     z    = A'*y; nprodAt=nprodAt+1;   % z = A'y throughout
  end

  % Some inactive constraints may be infeasible. Set their state.
  [sL,sU] = infeasibilities(bl,bu,z);
  inactive = abs(state) ~= 1;
  state(inactive  &  sL > featol) = -2;
  state(inactive  &  sU > featol) = +2;

  % Because active variables are feasible by definition (they're on their
  % respective bnds, we only need to look at inactive vars to determine if
  % the current point is infeasible.
  infeasible = any(abs(state)==2);
  
  if infeasible && ~elastic
    phase1  = true;  % Will turn off elastic mode as soon as
    elastic = true;  % problem is feasible.
    if loglevel > 0
      fprintf(fid,...
         [' Infeasible warm start. Enter phase 1. ' ...
	       ' No. of infeasible constraints: %i\n' ] ,...
         sum(abs(state)==2));
    end
  end     

%-----------------------------------------------------------------------
% Main loop.
%-----------------------------------------------------------------------
  while 1

    [sL,sU,maxInf] = infeasibilities(bl,bu,z);
    infeasible = any(abs(state)==2);
    newInfeas = infeasible && (nBrks > 0 || newInfeas);

    %------------------------------------------------------------------
    % Compute objective gradient and multipliers.
    %------------------------------------------------------------------
    g = b - lambda*y;        % Steepest-descent direction

    if elastic
        if infeasible
            if newInfeas     % Need to recompute Ae
                e = sign(state);
                e(active) = 0;  % Exclude indices in working set
                Ae = A*e;
                nprodA = nprodA+1;
            end
            g = g - rho*Ae;  % Add in contribution from penalty term
        elseif phase1
            % BPdual was in "phase1" mode. The current iterate is now
            % feasible, so change back to in-elastic mode.
            assert(maxInf < featol);
            elastic = false;
            if loglevel > 0
                fprintf(fid,...
                    [' == infeas %7.1e < %7.1e featol.'...
                     ' Exit phase 1. == \n'], maxInf,featol);
            end
        end
    end

    %------------------------------------------------------------------
    % Compute multipliers x, search directions (dx,dy), and residual r.
    %------------------------------------------------------------------
    if isempty(R)
        condS = 1;
    else
        rmin  = min(diag(R));
        rmax  = max(diag(R));
        condS = (rmax / rmin);
    end

    if condS > 1e+10
        eFlag = EXIT_SINGULAR_LS;
        npad  = size(S,2) - size(x,1);
        % Pad x with enough zeros to make it compatible with S.
        x = [x; zeros(npad,1)];
    else
        [dx,dy] = newtonstep(S, R, g, x,lambda);
        x = x + dx;
    end
    r = b - S*x;

    %------------------------------------------------------------------
    % Print to log.
    %------------------------------------------------------------------
    yNorm = norm(y,2);
    rNorm = norm(r,2);
    xNorm = norm(x,1);
    [pObj,dObj,rGap] = objectives(x,y,active,b,bl,bu,lambda,rNorm,yNorm);
    nact  = length(active);                  % size of the working set

    if loglevel > 0
        if     q, svar = sprintf('%8i%1s', q, '-');
        elseif p, svar = sprintf('%8i%1s', p ,' ');
        else      svar = sprintf('%8s%1s',' ',' ');
        end
        if mod(itn,50)==0
           fprintf(fid,'\n');
           fprintf(fid,logH,'Itn','step','Add/Drop','Active','rNorm2',...
               'xNorm1','Objective','RelGap','Infeas','nBrk','condS');
           if homotopy
              fprintf(fid,' %9s','lambda');
           end
           fprintf(fid,'\n');
        end
        fprintf(fid,logB,itn,step,svar,nact,rNorm,xNorm,...
                dObj,rGap,maxInf,nBrks,condS);
        if homotopy
           fprintf(fid,' %9.3e',lambda);
        end
        fprintf(fid,'\n');
    end

    %------------------------------------------------------------------
    % Check exit conditions.
    %------------------------------------------------------------------
    if ~eFlag && homotopy && lambda <= lamFinal
        eFlag = EXIT_OPTIMAL;
    end
    if eFlag == EXIT_OPTIMAL
        if elastic && infeasible
            eFlag = false;
            rho = 2*rho;
            if rho > rhomax
                eFlag = EXIT_INFEASIBLE;
            end
            if loglevel > 0
                fprintf(fid,' increasing rho to %8.0e\n',rho);
            end
        end
    end
    if ~eFlag && ~elastic && infeasible
        eFlag = EXIT_INFEASIBLE;
    end
    if ~eFlag && rGap < gaptol && nact > 0
        eFlag = EXIT_SMALL_DGAP;
    end
    if ~eFlag && itn >= itnMax
        eFlag = EXIT_TOO_MANY_ITNS;
    end
    if ~eFlag && nact >= actmax
        eFlag = EXIT_ACTMAX;
    end

    % Call the user's callback function
    if callback
        reqexit = ...
            opts.callback(itn,x,active,str2double(svar),lambda,r,y,pObj,dObj);
        if reqexit
           eFlag = EXIT_REQUESTED;
        end
    end

    % If this is an optimal solution, trim multipliers before exiting.
    if eFlag == EXIT_OPTIMAL || eFlag == EXIT_SMALL_DGAP
        if trim == 1
            % Optimal trimming. lambdain may be different from lambda.
            % Recompute gradient just in case.
            g = b - lambdain*y;
            [x,S,R,active,state] = trimx(x,S,R,active,state,g,b,lambda,opts);
            numtrim = nact - length(active);
            nact    = length(active);
        elseif trim == 2
            % Threshold trimming.
            % Not yet implemented.
        end
    end

    % Act on any live exit conditions.
    if eFlag
        break
    end

    %-------------------------------------------------------------------
    % New iteration starts here.
    %-------------------------------------------------------------------
    itn = itn + 1;
    p = []; q = [];

    %-------------------------------------------------------------------
    % Compute search direction dz.
    %-------------------------------------------------------------------
    if homotopy
        [x,dy,dz,step,lambda,p] = ...
            htpynewlam(active,state,A,R,S,x,y,sL,sU,lambda,lamFinal);
        nprodAt = nprodAt+1;        % htpynewlam incurs one mat-vec product
    else
        if norm(dy,inf) < eps       % save a mat-vec prod if dy is 0
            dz = zeros(n,1);
        else
            dz = A'*dy;  nprodAt = nprodAt+1;
        end
    end

    %---------------------------------------------
    % Find step to the nearest inactive constraint
    %---------------------------------------------
    if elastic
%       angles = (dz ./ cNorms) / norm(dy);
%       angles = dz;
%       gdy = -g'*dy;
        gdy = -lambda*norm(dy)^2;
        [step,state,nBrks,p,hitBnd] = ...
            l2maxline(gdy,dy,dz,z,bl,bu,lambda,rho,state); %,dz);
        if nBrks == 0, step = 1; end
    elseif homotopy
        if isempty(p)
            hitBnd = false;
        else
            hitBnd = true;
            if p < 0
                p = -p;
                state(p) = -1;
            else
                state(p) = +1;
            end
        end
    else
        stepL = inf; stepU = inf;
        blockL = find(dz<-pivtol & state==0); % This set moves down toward bl
        blockU = find(dz> pivtol & state==0); % This set moves up   toward bu
        if any(blockL)
            [stepL,pL] = min( (sL(blockL)-tietol)./dz(blockL) );
            pp         = blockL(pL);
            stepL      = sL(pp)/dz(pp);
        end
        if any(blockU)
            [stepU,pU] = min(-(sU(blockU)-tietol)./dz(blockU) );
            pp         = blockU(pU);
            stepU      = -sU(pp)/dz(pp);
        end
        stepL   = max(stepL,0);      % Don't allow a negative step.
        stepU   = max(stepU,0);      % Don't allow a negative step.
        step    = min(stepL,stepU);
        step    = min(step,1);
        hitBnd  = step < 1;
        if hitBnd                    % We bumped into a new constraint
            if step==stepL
                p = blockL(pL); state(p) = -1;
            else
                p = blockU(pU); state(p) =  1;
            end
        end
    end

    y = y + step*dy;                 % Update dual variables.
    if mod(itn,50)==0
       z = A'*y;                     % Occasionally set z = A'y directly.
       nprodAt = nprodAt+1; 
    else
       z = z + step*dz;              % Mostly, update z.
    end

    if hitBnd
       zerovec(p) = 1;               % Extract a = A(:,p)
       a          = A*zerovec;
       nprodA     = nprodA+1;
       zerovec(p) = 0;
       R          = QRaddcol(S,R,a); % Update R
       S          = [ S       a ];   % Expand S, active
       active     = [ active; p ];
       x          = [ x     ; 0 ];

    elseif elastic && nBrks > 0
       % An unconstrained step passed through some breakpoints.
       % Do another iteration on same active set.

    else
       %------------------------------------------------------
       % step = 1.  We moved to a minimizer for the current S.
       % See if we can delete an active constraint.
       %------------------------------------------------------
       drop = false;
       if any(active)
          if elastic
              dropl = state(active)==-1 & (x >  opttol | x < -rho-opttol);
              dropu = state(active)==+1 & (x < -opttol | x >  rho+opttol);
              zfrel = bl(active) + featol < z(active);
              zfreu = z(active) < bu(active) - featol;
              dropf = abs(x) > opttol & zfrel & zfreu;
              dropa = dropl | dropu | dropf;
          else
              dropl = state(active)==-1 & x >  opttol;
              dropu = state(active)==+1 & x < -opttol;
              dropa = dropl | dropu;
          end
          drop = any(dropa);
       end

       % Delete an active constraint.
       % Two cases:
       % 1. elastic mode: an active constraint may be infeasible. When
       % releasing that constraint, it's necessary to set its state
       % accordingly.
       % 2. inelastic mode: all released constraints are feasible.
       if drop
          [xmax,qa]  = max(abs(x.*dropa));
          q          = active(qa);
          newInfeas  = false;
          if ~elastic
              state(q) = 0;
          elseif dropf(qa)
              state(q) = 0;
          elseif state(q) == -1
              if x(qa) > 0
                  state(q) =  0;
              else
                  state(q) = -2;  newInfeas = true;
              end
          elseif state(q) == +1
              if x(qa) < 0
                  state(q) =  0;
              else
                  state(q) = +2;  newInfeas = true;
              end
          end
          S(:,   qa) = [];
          active(qa) = [];
          x(qa,:)    = [];             % keeps x a col vector even if empty
          R          = QRdelcol(R,qa);
       else
          eFlag = EXIT_OPTIMAL;
       end
    end % if step < 1
  end % while 1
%-----------------------------------------------------
% end main loop
%-----------------------------------------------------

% Gather exit data.
  inform.itns    = itn;
  inform.time    = toc(selftime);
  inform.stat    = eFlag - 1;
  inform.exitmsg = exit_msg{eFlag};
  inform.nprodA  = nprodA;
  inform.nprodAt = nprodAt;
  inform.rho     = rho;
  inform.rNorm   = rNorm;
  inform.xNorm   = xNorm;
  inform.opts    = opts;
  inform.numtrim = numtrim;

  if loglevel > 0
     fprintf(fid,'\n EXIT BPdual -- %s\n\n',exit_msg{eFlag});
     fprintf(fid,' %-20s: %8i %5s','No. significant nnz',sparsity(x),'');
     fprintf(fid,' %-20s: %8i\n','Products with A',nprodA);
     fprintf(fid,' %-20s: %8i %5s','No. trimmed nnz',numtrim,'');
     fprintf(fid,' %-20s: %8i\n','Products with At',nprodAt);
     fprintf(fid,' %-20s: %8.1e %5s','Solution time (sec)',inform.time,'');
     fprintf(fid,'\n');
  end

%---------------------------------------------------------------------
end % function BPdual.
%---------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRIVATE FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------------------------------------
% Compute a Newton step.  This is a step to a minimizer of the EQP
%
%   min   g'dy + 1/2 lambda dy'dy   subj to   S'dy = 0.
%
% The optimality conditions are given by
%
%   [ -lambda I   S ] [dy] = [ h ], with  h = b - lam y - S x, 
%   [     S'        ] [dx]   [ 0 ]
%
% where x is an estimate of the Lagrange multiplier. Thus, dx solves
% min ||S dx - h||. On input, g = b - lam y.
% ---------------------------------------------------------------------
function [dx,dy] = newtonstep(S, R, g, x, lambda)
[m,n] = size(S);
if m==0 || n==0
    dx = zeros(0,1);
    dy = g/lambda;          % Steepest descent
    return
end
h = g - S*x;
[dx,dr] = csne(R,S,h);      % LS problem  min ||S dx - h||.
if m > n                    % Overdetermined system.
    dy = dr/lambda;         % dy is the scaled residual.
else                        % System is square or underdetermined;
    dy = zeros(m,1);        % Anticipate that the residual is 0.
end

end % function newtonstep

%---------------------------------------------------------------------
% Compute the primal and dual objective values, and the duality gap:
%
%    DP:  minimize_y   - b'y  +  1/2 lambda y'y
%         subject to   bl <= A'y <= bu
%
%    PP:  minimize_x   bl'neg(x) + bu'pos(x) + 1/2 lambda y'y
%         subject to   Ax + lambda y = b.
%---------------------------------------------------------------------
function [pObj,dObj,rGap] = objectives(x,y,active,b,bl,bu,lambda,rNorm,yNorm)

bigNum = 1e20;

if isempty(x)
    blx = 0;
    bux = 0;
else
    blx = bl(active);
    blx(blx < -bigNum) = 0;
    blx = blx' * min(x,0);
    
    bux = bu(active);
    bux(bux >  bigNum) = 0;
    bux = bux' * max(x,0);
end

if lambda > eps
    pObj = blx + bux + lambda\rNorm^2/2;   % primal objective
else
    pObj = blx + bux;
end
dObj  = lambda*yNorm^2/2 - b'*y;           % dual   objective

maxpd = max(max(1,pObj),dObj);
rGap  = abs(pObj+dObj)/maxpd;              % relative duality gap

end % function objectives

%---------------------------------------------------------------------

function [sL,sU,maxInf] = infeasibilities(bl,bu,z)
%INFEASIBILITIES  Compute the infeasibility of z relative to bounds bl/bu.
%
%  (bl - z) <= 0  implies  z is   feasible wrt to lower bound
%           >  0  ...           infeasible ...
%  (z - bu) <= 0  implies  z is   feasible wrt to lower bound.
%           >  0  ...           infeasible ...

sL = bl - z;
sU = z - bu;
if nargout > 2
   maxInf = max(0,max(max(sL),max(sU)));
end
end % function infeasibilities
