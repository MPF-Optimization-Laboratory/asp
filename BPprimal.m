function [x, y, inform] = BPprimal( A, b, delta, lambda, x, options, inform )

%        [x, y, inform] = BPprimal( A, b, delta, lambda, x, options, inform );
%
% Solve the regularized basis pursuit denoise (BPDN) problem
%
%    minimize    ||x||_1  +  1/2 delta x'x  +  1/2 lambda y'y
%       x,y
%    subject to   Ax + lambda y = b,
%
% where
% A          is an m x n linear operator (explicit or implicit).
% b          is an m-vector.
% delta      is a positive scalar (typically small, say 1e-4).
% lambda     is a positive scalar.
% x          is an estimate of the solution.
% options    is a structure of options from l1Set. Any unset options
%            are set to their default value; set options=[] to use all
%            default values.
% inform     is the output from a previous call to BPprimal. Use this
%            argument to warm-start the algorithm; or set inform=[]
%            for a cold start.
%
% OUTPUT ARGUMENTS:
% x          is the solution.
% y          is a scaled residual.  The true residual is r = delta y = b - Ax.
% inform     is a structure with exit information:
% .stat      Exit conditions:
%              0 = optimal,
%              1 = too many major its,
%              2 = too many minor its.
% .mjrIts    Number of major iterations
% .mnrIts    Number of minor iterations
% .time      cpu time used
%
% AUTHORS
% Michael P. Friedlander               Michael A. Saunders
% University of British Columbia       Stanford University
% mpf@cs.ubc.ca                        saunders@stanford.edu
%
% 09 Nov 2007: Renamed from lsl1_single_lambda to BPprimal.
% 09 Aug 2008: delta and y introduced into problem statement
%              so that both delta and lambda may be regarded
%              as perturbations to the original BP problem.
% 12 Aug 2008: Implemented BFGS approximation to reduced Hessian:
%                 R'R ~= S'S + delta*lambda*I.
%
% $Id: BPprimal.m 515 2010-04-20 16:48:10Z mpf $
%----------------------------------------------------------------------
REVISION = '$Revision: 515 $';
DATE     = '$Date: 2010-04-20 09:48:10 -0700 (Tue, 20 Apr 2010) $';
REVISION = REVISION(11:end-1);
DATE     = DATE(8:26);

%%---------------------------------------------------------------------
%% Set default input arguments.
%%---------------------------------------------------------------------
  explicitA = ~isa(A,'function_handle');

  if nargin < 2
     error('At least two arguments are required');
  end
% If A is not explicit, then x must be given.
  nogivenx = nargin < 5  || isempty(x); 
  if explicitA
     n = size(A,2);
  elseif nogivenx
     error('A is implicit, so x must be given.');
  else
     n = length(x);
  end
  if nargin < 6, options = []; end
  if nargin < 7, inform  = []; end
  m = length(b);
  
%----------------------------------------------------------------------
% Grab input options and set defaults where needed. 
%----------------------------------------------------------------------
  if isempty(options)
     options = as_setparms;
  end
  fid          = options.fid;
  maxMjr       = options.majiters*max(n,m);
  maxMnr       = options.miniters;
  mprice       = options.mprice;
  feaTol       = options.featol;
  optTol       = options.opttol;
  subTol       = options.subtol;
  lbfgsmem     = options.lbfgsmem;
  callback     = isa(options.callback,'function_handle');
  
  LSmethod     = options.LSmethod;
  LS.BFGS      = strcmpi(LSmethod, 'BFGS'    );
  LS.LBFGS     = strcmpi(LSmethod, 'LBFGS'   ); % Limited-memory BFGS
  LS.BFGSdiag  = strcmpi(LSmethod, 'BFGSdiag');
  LS.LSQR      = strcmpi(LSmethod, 'LSQR'    );
  LS.PCG       = strcmpi(LSmethod, 'PCG'     );
  LS.QlessQR   = strcmpi(LSmethod, 'QlessQR' );
  LS.QR        = strcmpi(LSmethod, 'QR'      );
  LS.steepest  = strcmpi(LSmethod, 'Steepest'); % Steepest descent
  LS.nonlinCG  = strcmpi(LSmethod, 'nonlinCG'); % Nonlinear CG
  if LS.QR
     if ~isnumeric(A)
        error('QR requires A to be a matrix!');
     end
  elseif LS.QlessQR
     zerovec = zeros(n,1); % Needed to extract columns from A
  end

%----------------------------------------------------------------------
% Initialize various variables.
%----------------------------------------------------------------------
  tic;                    % Start your watches!
  itn       = 0;          % Major iteration counter
  minIts    = 0;          % Minor iteration counter
  nLine     = 0;          % No. of linesearch iterations
  sTol      = inf;        % Forces a price on the very first iteration
  condH     = 0;          % Lower bound on condition of approx reduced H
  condHmax  = 1e+8;       % Limit on condH before reset
  step      = 0;          % Step length
  newS      = false;      % Indicates if a new superbasic entered
  nprodA    = 0;          % Count products A*x
  nprodAt   = 0;          % Count products A'*y
  nSdel     = 0;          % No. deleted from S by knock-out
  beta      = sqrt(delta*lambda);      % Damping parameter for LS problem
  dscale    = sqrt(1/lambda + delta);  % Diagonal scaling for R,
                                       % where R'R ~= reduced Hessian
  xdual     = zeros(n,1); % Used to compute duality gap
  
% Exit flags.
  EXIT_OPTIMAL         = 1;
  EXIT_TOO_MANY_MAJORS = 2;
  EXIT_TOO_MANY_MINORS = 3;
  EXIT_LINESEARCH_FAIL = 4;
  EXIT_ZERO_STEP       = 5;
  EXIT_REQUESTED       = 6;
  exit_msg = {
      'Optimal'
      'Too many major iterations'
      'Too many minor iterations'
      'Linesearch failed'
      'Zero step'
      'Exit requested'
             };
  
%-----------------------------------------------------------------------
% Initialize cold/warm/hot-start variables.
%-----------------------------------------------------------------------
  coldstart = nogivenx;
  warmstart = ~coldstart && isempty(inform);
  S = [];                       % Storage for superbasic columns of A
  R = [];                       % BFGS or QR triangle
  D = [];                       % BFGS diagonal
  H = [];                       % L-BFGS structure
  if coldstart
  % --------------------------------------------------------------------
  % Coldstart. No superbasics, no Hessian information.
  % --------------------------------------------------------------------
     iS = zeros(1,0);
     nS = 0;
     
  elseif warmstart
  % --------------------------------------------------------------------
  % Warmstart. Infer active set from x.  No Hessian information.
  % --------------------------------------------------------------------
     iS = find(abs(x) > feaTol)';
     nS = length(iS);

     if LS.QlessQR
     %  QlessQR: Assemble superbasic columns and corresponding QR.
        S = zeros(n,nS);
        for j=iS
            zerovec(j) = 1;
            Snew   = Aprod(zerovec,1);         % Grab j-th column from A
            R      = QRaddcol(Snew,R,a,beta);  % Update R
            S(:,j) = Snew;                     % Add new column to S
            zerovec(j) = 0;
        end
        nprodA = nprodA + nS;
        
     elseif LS.BFGS
        R = dscale*eye(nS);
        
     elseif LS.BFGSdiag
        D = repmat(dscale,nS,1);
     end
     
  else
  % --------------------------------------------------------------------
  % Hotstart. Grab active-set and Hessian info from inform structure. 
  % --------------------------------------------------------------------
     iS   = inform.iS;
     nS   = length(iS);
     sTol = inform.sTol;
     if LS.QlessQR
        error('Hotstart for QlessQR not yet implemented');
        
     elseif LS.BFGS
        if ~isempty(inform.R)
            R = inform.R;
        elseif ~isempty(inform.D)
            R = diag(inform.D);
        end
        
     elseif LS.BFGSdiag
        if ~isempty(inform.R)
           D = diag(inform.R);
        elseif ~isempty(inform.D)
           D = inform.D;
        else
           D = repmat(dscale,nS,1);
        end
     end
  end
  iN      = (1:n);             % iN is the complement of iS
  iN(iS)  = [];

%---------------------------------------------------------------------
% Initialize method-dependent variables.
%---------------------------------------------------------------------
  atol1       = 1e-5;   % For lsqr 
  atolmin     = 1e-15;
  atolcurrent = atol1;
  itnCG       = 0;
  printCG     = false;
  if LS.LBFGS
     H = lbfgsinit(nS,lbfgsmem,1/dscale);
  end
  
%-----------------------------------------------------------------------
% Initialize Aty, gZ, etc.
%-----------------------------------------------------------------------
  if coldstart
     x      = zeros(n,1);
     y      = b/lambda;
  else
     Ax     = Aprod(x,1);
     nprodA = nprodA + 1;
     y      = (b - Ax)/lambda;
  end

  Aty       = Aprod(y,2);
  nprodAt   = nprodAt + 1;

  if nS==0
     gZ = zeros(nS,1);
  else
     xS = x(iS);
     gS = sign(xS) + delta*xS;
     gZ = gS - Aty(iS);
  end
  gZ0 = gZ;
  dxS = zeros(nS,1);
  dx  = zeros(n ,1);
  
%-----------------------------------------------------------------------
% Iteration log.
%-----------------------------------------------------------------------
  logHead = sprintf('\n %4s %6s %6s %8s %7s %8s %12s %7s %7s %7s %5s %5s %5s',...
                    'Itn','nS','Sadd','Optimal','gZNorm','Rel Gap','Primal Obj',...
                    'Resid','xNorm1','step','nLine','nSdel','condH');
  logBody1 = '\n %4i %6i %6i %8.1e %7.1e %8.1e %12.6e %7.1e %7.1e %7.1e %5i %5i %5.0e';
  logBody2 = '\n %4i %6i %6s %8s %7.1e %8.1e %12.6e %7.1e %7.1e %7.1e %5i %5i %5.0e';
  fprintf(fid,'\n');
  fprintf(fid,' =================================================\n');
  fprintf(fid,' %s  v.%s (%s)\n', mfilename, REVISION, DATE);
  fprintf(fid,' =================================================\n');
  fprintf(fid,'\n');
  fprintf(fid,' Parameters:');
  fprintf(fid,'\n %20s = %7i'  , 'No. of equations   ' , m       );
  fprintf(fid,'   %20s = %7.1e', 'Feasibility tol    ' , feaTol  );
  fprintf(fid,'\n %20s = %7i'  , 'No. of variables   ' , n       );
  fprintf(fid,'   %20s = %7.1e', 'Optimality  tol    ' , optTol  );
  fprintf(fid,'\n %20s = %7i'  , 'Max no of major its' , maxMjr  );
  fprintf(fid,'   %20s = %7.1e', 'Subspace    tol    ' , subTol  );
  fprintf(fid,'\n %20s = %7.1e', 'delta              ' , delta   );
  fprintf(fid,'   %20s = %7i'  , 'Multiple price     ' , mprice  );
  fprintf(fid,'\n %20s = %7.1e', 'lambda             ' , lambda  );
  fprintf(fid,'   %20s = %7s'  , 'LS method          ' , LSmethod);
  fprintf(fid,'\n %20s   %7s'  , '                   ' , ''      );
  if LS.LBFGS
  fprintf(fid,'   %20s = %7i'  , 'L-BFGS vectors     ' , lbfgsmem);
  end
  fprintf(fid,'\n');
  fprintf(fid,logHead);
  if LS.LSQR
     fprintf(fid,' %6s', 'LSQR')
  elseif LS.PCG
     fprintf(fid,' %6s', 'PCG' )
  end
  
%-----------------------------------------------------------------------
% Major iterations
%-----------------------------------------------------------------------
  while 1
    
  %---------------------------------------------------------------------
  % Compute objective, optimality, etc.
  %---------------------------------------------------------------------
    if itn > 0  &&  mod(itn, 10000)==0
       Ax     = Aprod(x,1);
       nprodA = nprodA + 1;
       y      = (b - Ax)/lambda;
    end

    xNorm1 = norm(x,1);   % Could be xS here, except some xN may be nonzero
    xNorm2 = norm(x,2);   % Ditto
    yNorm2 = norm(y,2);
    yNormI = norm(y,inf);
    rNorm2 = lambda*yNorm2;
    gZNorm = norm(gZ,inf);

    % Compute xdual, the largest dual feasible x for current y.
    J         = Aty >= 0;
    xdual( J) = max( Aty( J)-1, 0 );
    xdual(~J) = min( Aty(~J)+1, 0 );
    xdual     = xdual/delta;
    xdNorm2   = norm(xdual,2);

    % Compute objectives and duality gap for
    % primal feasible (x,y) and dual feasible (xdual,y).
    pobj   = xNorm1 + 0.5*(delta*xNorm2^2  + lambda*yNorm2^2);
    dobj   = b'*y   - 0.5*(delta*xdNorm2^2 + lambda*yNorm2^2);
    rgap   = (pobj - dobj) / (1+abs(pobj)+abs(dobj));
    
  %---------------------------------------------------------------------
  % Print log.
  %---------------------------------------------------------------------
    if newS
       fprintf(fid,logBody1,itn,...
               nS,iSadd(1),rgmax,gZNorm,rgap,pobj,rNorm2,xNorm1,step,...
               nLine,nSdel,condH);
    else
       fprintf(fid,logBody2,itn,...
               nS,'     ','    ',gZNorm,rgap,pobj,rNorm2,xNorm1,step,...
               nLine,nSdel,condH);
    end

    if printCG, fprintf(fid,' %6i', itnCG); end
    if LS.LSQR
       if itnCG==1
          if atolcurrent > atolmin+eps
             atolcurrent = max( 0.1*atolcurrent, atolmin );
             fprintf(fid,'\n atol reduced to %8.1e', atolcurrent)
          end
       end
    end
    
  %---------------------------------------------------------------------
  % Check exit conditions.
  %---------------------------------------------------------------------
    % Call the user's callback function
    if callback
        reqexit = options.callback(itn,x,rNorm2);
        if reqexit
           flag = EXIT_REQUESTED;
        end
    elseif rgap <= optTol
       flag = EXIT_OPTIMAL;
       break
    elseif itn >= maxMjr
       flag = EXIT_TOO_MANY_MAJORS;
       break
    elseif itnCG >= maxMnr
       flag = EXIT_TOO_MANY_MINORS;
       break
    end
    
  %---------------------------------------------------------------------
  % New major iteration starts here.
  %---------------------------------------------------------------------
    itn     = itn + 1;
    newS    = false;
    nSdel   = 0;
    gZadd   = zeros(0,1);
    nSadd   = 0;          rgmax   = [];     % if we don't call BPprice.
    updateR = false;      printCG = false;
    updateD = false;      itnCG   = 0;
    updateH = false;
    
  %---------------------------------------------------------------------
  % Price (only if we're near optimal on the current superbasics).
  %---------------------------------------------------------------------
  % if nS >= m, sTol = optTol; end % With delta>0, this no longer seems relevant

    if gZNorm < sTol % && nS <= m

       gTol  = optTol*(1 + yNormI);            % Relative optimality tolerance

       [gZadd,iNdel,rgmax] = BPprice( Aty, x, delta, iN, mprice, gTol, subTol );

       iSadd = iN(iNdel);
       nSadd = length(iSadd);
       sTol  = max( subTol*rgmax, gTol ); % Next target gZNorm

       %----------------------------------------------------------------
       % Knock-out.  Move small variables from S to the end of N.
       %----------------------------------------------------------------
 %     Sdel  = find( abs(x(iS)) < feaTol  &  abs(gZ) < sTol )';
       Sdel  = find( abs(x(iS)) < feaTol  &  abs(gZ) < gTol )';
       nSdel = length(Sdel);
       if nSdel
          x(iS(Sdel)) = 0;
          dxS(Sdel)   = [];
          [iN,iS,nS,gZ,gZ0,R,D,H,S,gZNorm] = knockout(LS,iN,iS,Sdel,gZ,gZ0,R,D,H,S);
       end

       if rgmax <= gTol
          if gZNorm  <= gTol
             flag = EXIT_OPTIMAL;
             if nSdel
                fprintf(fid,'\n Final knock-out: %i', nSdel)
             end
             break
          end
       elseif nSadd
          newS      = true;
          dxS       = [dxS; zeros(nSadd,1)];
          iS        = [iS         iSadd   ];
          gZ0       = [gZ0;       gZadd   ];
          gZ        = [gZ ;       gZadd   ];
          gZNorm    = max(gZNorm,rgmax);
          iN(iNdel) = [];
          nR        = nS;
          nS        = length(iS);

          if LS.BFGS
             R = [ R                  zeros(nR,nSadd)
                   zeros(nSadd,nR)  dscale*eye(nSadd) ];
               
          elseif LS.BFGSdiag
             D = [         D
                   repmat(dscale,nSadd,1) ];
               
          elseif LS.LBFGS
             H = lbfgsadd(H,nSadd,1/dscale);
             
          elseif LS.QlessQR
             for j = iSadd
                 zerovec(j) = 1;
                 Snew = Aprod(zerovec,1);         % Grab j-th column from A
                 R    = QRaddcol(S,R,Snew,beta);  % Update R
                 S    = [S Snew];                 % Add new column to S
                 zerovec(j) = 0;
             end
             nprodA = nprodA + nSadd;
             
          end
       end
    end
    xS  = x(iS);
 
  %---------------------------------------------------------------------
  % Compute search directions dx, dy.
  % Z'HZ dxS = -gZ,  where Z'HZ ~= R'R (BFGS) or D (BFGSdiag).
  % This is equivalent to the LS problem
  %   min ||[S; beta*I]w - [y; -gS/beta]||,  dxS = lambda*w.
  %
  % Then  dy = -1/lambda dxS.
  %---------------------------------------------------------------------
    if LS.steepest || (newS && ~(LS.QR || LS.QlessQR))  % Try steepest descent after Price
      dxS = -gZ;

    elseif LS.BFGS
      v   = R'\ gZ;
      dxS = R \ -v;
      updateR = true;

    elseif LS.BFGSdiag
      dxS = -gZ ./ D;
      updateD = true;
      
    elseif LS.LBFGS
      dxS = lbfgshprod(H,-gZ);
      updateH = true;
      
    elseif LS.LSQR  % Use LSQR to solve the same problem as for LS.QR
      sgnxS = sign(xS);
      if newS
         sgnxS(end-nSadd+1:end) = -sign(gZadd);
      end
      gS         = sgnxS + delta*xS;
      s          = [   y
                    -gS/beta];
      lsqrSIprod = @(x,mode) SIprod( x,mode,@Aprod,m,n,iS,beta );
      damp       = 0;
      atol       = atolcurrent;
      atol       = min( atol, 1e-3 );
      btol       = atol;
      conlim     = 1e+6;
      itnlim     = 10*nS;    % Could be related to maxMnr
      show       = 0;

      [ w, istop, itnCG, r1norm, r2norm, SInorm, SIcond ]...
         = lsqr( m+nS, nS, lsqrSIprod, s, ...
                 damp, atol, btol, conlim, itnlim, show );

      dxS        = lambda*w;
      nprodA     = nprodA  + itnCG;
      nprodAt    = nprodAt + itnCG + 1;
      minIts     = minIts  + itnCG;
      condH      = SIcond;           % For printing
      printCG    = true;
      
    elseif LS.PCG
      pcgStSIprod = @(x) StSIprod( x,@Aprod,n,iS,beta );
      M1     = []; M2 = [];
      pcgtol = min(.5,sqrt(gZNorm));
      itnlim = 10*nS;

      % Don't use last search direction as starting point after Price.
      if newS
         [ w, flagPCG, relresPCG, itnCG ] = ...
             pcg( pcgStSIprod, -gZ, pcgtol, itnlim, M1, M2 );
      else          
         w = dxS / lambda;      
         [ w, flagPCG, relresPCG, itnCG ] = ...
             pcg( pcgStSIprod, -gZ, pcgtol, itnlim, M1, M2, w );
      end
      
      dxS     = lambda*w;
      nprodA  = nprodA  + itnCG;
      nprodAt = nprodAt + itnCG;
      minIts  = minIts  + itnCG;
      printCG = true;
      
    elseif LS.QlessQR
      w   = R \ (R'\-gZ);
      dxS = lambda*w;
      
    elseif LS.QR
      sgnxS = sign(xS);
      if newS
         sgnxS(end-nSadd+1:end) = -sign(gZadd);
      end
      gS  = sgnxS + delta*xS;
      bI  = beta*speye(nS);
      if ~issparse(A), bI = full(bI); end
      SI  = [A(:,iS)
               bI   ];
      s   = [   y
             -gS/beta];
      w   = SI\s;    % LS problem (dense or sparse)
      dxS = lambda*w;
      
    elseif LS.nonlinCG
      beta = (gZ'*gZ - gZ'*gZ0) / (gZ0'*gZ0);
      beta = max(beta,0);
      if beta > 100; beta = 0; end
      dxS  = -gZ + beta*dxS;
    end

    if LS.QlessQR
       dy = -S*dxS/lambda;
    else
       dx(iN) = 0;
       dx(iS) = dxS;
       dy     = - Aprod(dx,1)/lambda;  % B*dy = - S*dx, where B = lambda*I
       nprodA = nprodA + 1;
    end
 
  %---------------------------------------------------------------------
  % Linesearch.
  %---------------------------------------------------------------------
    [step,nLine,LNstat,Sdel,Ssgn] = l2l1line(y,dy,xS,dxS,lambda,delta,1);

    if LNstat == 4  % unbounded step
       flag = EXIT_LINESEARCH_FAIL;
       break
    end
    
    if step==0
       if LS.LSQR
          if atolcurrent > 1.1*atolmin
             atolcurrent = max( 0.1*atolcurrent, atolmin );
             fprintf(fid,'\n atol reduced to %8.1e', atolcurrent)
             continue
          elseif nS <= 1000
             LS.LSQR = false;
             LS.QR   = true;
             fprintf(fid,'\n Switching to QR')
             continue
          end
       end
       
       flag = EXIT_ZERO_STEP;
       break
    end
    
    xS    = xS + step*dxS;
    x(iS) = xS;
    y     = y  + step*dy;

  %---------------------------------------------------------------------
  % Update Aty, gZ, R or D.
  %---------------------------------------------------------------------
    Aty     = Aprod(y,2);
    nprodAt = nprodAt + 1;

    gZ0     = gZ;
    gS      = sign(xS);
    if LNstat == 1
       gS(Sdel) = Ssgn;
    end
    gS      = gS + delta*xS; 
    gZ      = gS - Aty(iS);

  % Update quasi-Newton reduced Hessian
    if LS.BFGS && updateR
       [R,noupdate] = bfgsupdate(R, step, dxS, gZ0, gZ, v);
       if noupdate
          fprintf(fid,'\n  W: bfgs update skipped');
       end
       diagR  = abs(diag(R));
       condH  = ( max(diagR) / min(diagR) )^2;
       if condH > condHmax
          R  = dscale*eye(nS);
          fprintf(fid,'\n  W: Reset R')
       end
    elseif LS.BFGSdiag && updateD
       D     = bfgsdiag  (D, step, dxS, gZ0, gZ, dscale);
       condH = max(D) / min(D);   % Diagonals of D are positive
       if condH > condHmax
          D  = repmat(dscale,nS,1);
          fprintf(fid,'\n  W: Reset D')
       end
    elseif LS.LBFGS && updateH
        [H,noupdate] = lbfgsupdate(H, step, dxS, gZ0, gZ);
        if noupdate
           fprintf(fid,'\n  W: l-bfgs update skipped');
        end           
    end
    
  % The search stopped at a breakpoint. Delete corresponding superbasic.
    if LNstat == 1   
       gZ(Sdel) = 0;  % Force the next line to include this single Sdel
       Sdel     = find( abs(xS) < feaTol  &  abs(gZ) < sTol )';
       if any(Sdel)
          nSdel       = nSdel + length(Sdel);
          x(iS(Sdel)) = 0;
          dxS(Sdel)   = [];
          [iN,iS,nS,gZ,gZ0,R,D,H,S,gZNorm] = knockout(LS,iN,iS,Sdel,gZ,gZ0,R,D,H,S);
       end
    end

  end % while ~err

%---------------------------------------------------------------------
% Gather exit data.
%---------------------------------------------------------------------
  inform.majIts  = itn;
  inform.minIts  = minIts;
  inform.time    = toc;
  inform.stat    = flag;
  inform.iS      = iS;
  inform.R       = R;
  inform.D       = D;
  inform.sTol    = sTol;
  inform.y       = y;
  inform.nprodA  = nprodA;
  inform.nprodAt = nprodAt;
  
  if fid
  fprintf(fid,'\n\n EXIT BPprimal -- %s\n',exit_msg{flag});
  fprintf(fid,' %-20s: %8i %5s','No. significant nnz',sparsity(x),'');
  fprintf(fid,' %-20s: %8i\n','Products with A',nprodA);
  fprintf(fid,' %-20s: %8.1e %5s','Solution time (sec)',inform.time,'');
  fprintf(fid,' %-20s: %8i\n','Products with At',nprodAt);
  fprintf(fid,'\n');
  end
  
  return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED FUNCTIONS.  These share some vars with workspace above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = Aprod(x,mode)
   if mode == 1
      if   explicitA, z = A*x;
      else            z = A(x,1);
      end
   elseif mode == 2
      if   explicitA, z = A'*x;
      else            z = A(x,2);
      end
   else
      error('Wrong mode!');
   end
end % function Aprod

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % function BPprimal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRIVATE FUNCTIONS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [iN,iS,nS,gZ,gZ0,R,D,H,S,gZNorm] = knockout(LS,iN,iS,Sdel,gZ,gZ0,R,D,H,S)

% Important: we assume that Sdel is sorted (ascending) on input.
% 03 Sep 08: Found a bug: was deleting columns in R using
%            fliplr(find(Sdel)) rather than fliplr(Sdel).

  Nnew     = iS(Sdel);
  iN       = [iN Nnew];
  iS(Sdel) = [];
 gZ0(Sdel) = [];
  gZ(Sdel) = [];
  gZNorm   = norm(gZ,inf);  % Note: gZNorm may become small

  if LS.BFGS
    for j = fliplr(Sdel)
      R = QRdelcol(R,j);
    end
    
  elseif LS.QlessQR
    for j = fliplr(Sdel)
      R = QRdelcol(R,j);
    end     
    S(:,Sdel) = [];
    
  elseif LS.BFGSdiag
    D(Sdel) = [];
    
  elseif LS.LBFGS
    H = lbfgsdel(H,Sdel);
    
  end
  nS  = length(iS);
end % function knockout

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = SIprod( x, mode, Aprod, m, n, iS, beta )

% Used by lsqr to solve w = SI\s;
% Implements products for the matrix SI = [ S; beta*I ], where S = A(:,iS).
  if mode==1   % y = [S; beta I]*x
    xn     = zeros(n,1);
    xn(iS) = x;
    y      = [ Aprod(xn,1)
                beta*x    ];
  else         % y = S'*x
    nS    = length(iS);
    yn    = Aprod(x(1:m),2);
    y     = yn(iS) + beta*x(m+1:m+nS);
  end
end % function SIprod

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = StSIprod( x, Aprod, n, iS, beta )

% Used by pcg to solve  [S'S + beta^2I] w = -gZ
    
  xn     = zeros(n,1);
  xn(iS) = x;
  w      = Aprod(xn,1);        % w =   S x
  z      = Aprod( w,2);        % z = S'S x
  y      = z(iS) + beta^2*x;   % y = S'S x + beta^2 x
    
end % function StSIprod

