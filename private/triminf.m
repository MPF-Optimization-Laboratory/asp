function [active,state,S,R] = triminf(active,state,S,R,bl,bu,g,opts)
%TRIMINF   Trim working constraints with "infinite" bounds.

logB = '%4i %4i %4i %10.3e %10.e %10.e %10.e\n';
bigbnd = 1e10;
nact = length(active);

% Generate a list of constraints to trim.
tlistbl = find( state(active) == -1 & bl(active) < -bigbnd );
tlistbu = find( state(active) == +1 & bu(active) > +bigbnd );
tlist   = [tlistbl; tlistbu];

if isempty(tlist)
   return
end

for k = 1:length(tlist)
    q  = tlist(k);
    qa = active(q);         % Index of active constraint to delete
    nact = nact - 1;

    S(:   ,qa) = [];        % Delete column from S
    active(qa) = [];        % Delete index from active set
    R = QRdelcol(R,qa);     % Recompute new QR factorization
    state(q) = 0;           % Mark constraint as free
    x = csne(R,S,g);        % Recompute multipliers

    rNorm = norm(b - S*x);
    xNorm = norm(x,1);
    printf(opts,logB,k,q,nact,bl(q),bu(q),rNorm,xNorm);
end
    

