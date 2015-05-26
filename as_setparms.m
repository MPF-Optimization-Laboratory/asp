function defopts = as_setparms(opts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin default options.  Use only lower-case fields.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generic options
defopts.featol     = 5e-05  ;
defopts.gaptol     = 1e-06  ;
defopts.opttol     = 1e-05  ;
defopts.callback   =    []  ;
defopts.fid        =     1  ;
defopts.minfid     =     0  ;
defopts.loglevel   =     1  ;

% BPprimal options
defopts.subtol     =   0.2  ;
defopts.mprice     =    10  ;
defopts.lbfgsmem   =     5  ;
defopts.LSmethod   = 'bfgs' ;

% BPdual options
defopts.pivtol     = 1e-12  ;
defopts.homotopy   = false  ;
defopts.majiters   =    10  ;
defopts.miniters   = 10000  ;
defopts.actmax     =   inf  ;  % Max. no. of allowable active constraints
defopts.elastic    = false  ;
defopts.phase1     =  true  ;
defopts.rhomax     = 1e+06  ;
defopts.rho        = 1e+00  ;
defopts.trim       =     1  ;

% Reweighted BP options
defopts.rwfid      =     1  ;
defopts.rwtol      = 1e-01  ;
defopts.rwfeaTol   = 1e-06  ;
defopts.rwitns     =    25  ;
defopts.rwtopy     =  true  ;

% Sequential CS options
defopts.num_additional_measurements = 1;
defopts.tol_seqcs_convergence = 1e-3;

% OMP options
defopts.ompitns    =     1  ;   % iterations: k*m

% Generic parameters
defopts.lambdamin  = sqrt(eps); % minimum allowable lambda

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End default options.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Quick return if no user opts
if nargin == 0 || isempty(opts) 
    return
end

% List of valid field names
vfields = fieldnames( defopts );

% Grab valid fields from user's opts
for i = 1:length(vfields)
    field = vfields{i};
    if isfield( opts, field );
        defopts.(field) = opts.(field);
    end
end
