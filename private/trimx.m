function [x,S,R,active,state] = trimx(x,S,R,active,state,g,b,lambda,opts)
%TRIMX   Trim unneeded constraints from the active set.
%
% TRIMX assumes that the current active set is optimal, ie,
% 1. x has the correct sign pattern, and
% 2. z := A'y is feasible.
%
% Keep trimming the active set until one of these conditions is
% violated. Condition 2 isn't checked directly. Instead, we simply check if
% the resulting change in y, dy, is small.  It would be more appropriate to
% check that dz := A'dy is small, but we want to avoid incurring additional
% products with A. (Implicitly we're assuming that A has a small norm.)

logH = '\n %4s  %8s %8s  %7s  %10s  %10s  %10s\n';
logB = ' %4iT %8.1e %8i- %7i  %10.4e  %10.4e  %10.4e\n';
k = 0;
opttol = opts.opttol;
featol = opts.featol;
nact = length(active);
xabs = abs(x);
[xmin,qa] = min(xabs);
gNorm = norm(g,inf);

while xmin < opttol

    e = sign(x.*(xabs > opttol));  % Signs of significant multipliers
    q = active(qa); % Index of the corresponding constraint
    a = S(:,qa);    % Save the col from S in case we need to add it back. 
    xsmall = x(qa); % Value of candidate multiplier

    % Trim quantities related to the small multiplier.
    e(qa,:)    = []; % Ensure that e remains a column vector even if empty
    S(:   ,qa) = [];
    active(qa) = [];
    R = QRdelcol(R,qa);

    % Recompute the remaining multipliers and their signs.
    [x,dy] = csne(R,S,g);           % min ||g - Sx||_2
    xabs = abs(x);
    et = sign(x.*(xabs > opttol));
    dyNorm = norm(dy,inf)/lambda;   % dy = (g - Sx) / lambda
    
    % Check if the trimmed active set is still optimal
    if any( et ~= e ) || (dyNorm / max(1,gNorm) > featol)
       R      = QRaddcol(S,R,a);
       S      = [ S      a ]; %#ok<*AGROW>
       active = [active; q ];
       x      = csne(R,S,g);
       return
    end

    if opts.loglevel > 0 && mod(k,50)==0
       fprintf(opts.fid,logH,...
           'Itn','xSmall','Add/Drop','Active','rNorm2','xNorm1','dyNorm');
    end
    
    % The trimmed x is still optimal.
    k     = k + 1;
    nact  = nact - 1;
    state(q) = 0;               % Mark this constraint as free.
    rNorm = norm(b - S*x);
    xNorm = norm(x,1);
    [xmin,qa] = min(xabs);      % Grab the next candidate multiplier.
    if opts.loglevel > 0
        fprintf(opts.fid,logB,k,xsmall,q,nact,rNorm,xNorm,dyNorm);
    end
end % while
end % function trimx
