function inform = demo_seqcs
%DEMO_SEQCS  Demo the routine as_seqcs.
%
%   See also AS_SEQCS.

% BPdual Toolbox
% Copyright 2010, Michael P. Friedlander and Michael A. Saunders
% http://www.cs.ubc.ca/labs/scl/bpdual
%
% $Id: demo_seqcs.m 519 2010-04-21 16:36:30Z mpf $

% "New" way of resetting randn and rand
reset(RandStream.getDefaultStream);
   
% Problem dims
n  = 1024;
k  = 10;
    
% Sparse solution (x0) and signal (y)
p0 = randperm(n)'; p0 = p0(1:k);
x0 = zeros(n,1);
x0(p0) = sign(randn(k,1));
y  = dct(x0);

% Predetermine the order in which to sample the signal.
idx = randperm(n)';

% MB is a function that returns the composite k-by-n operator M*B
MB  = @(k)operator(@(w,mode)restricted_dct(n,idx(1:k),w,mode),k,n);

% My is a function that returns the sampled signal
My  = @(k)(y(idx(1:k)));

% Sequential CS for basis pursuit (ie, lambda = 0)
opts = as_setparms;
opts.num_additional_measurements = 10; % Grab 10 additional measurements at a time
lambda = 0;
[x,inform] = as_seqcs(MB,My,lambda,opts);

end % demo_seqcs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Private functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v = restricted_dct(n,idx,w,mode)
%RESTRICTED_DCT  Restricted DCT operator.
% If mode == 1, returns  v =  RB  *w;
% if mode == 2, returns  v = (RB)'*w,
% where R = eye(:,idx) and B is a DCT operator.
if mode == 1
   z = dct(w);
   v = z(idx);
elseif mode == 2
   v = zeros(n,1);
   v(idx) = w;
   v = idct(v);
end
end
