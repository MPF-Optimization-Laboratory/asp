function [dT,kT,rgmax] = BPprice( Aty, x, delta, iN, mprice, gTol, subTol )
%        [dT,kT,rgmax] = BPprice( Aty, x, delta, iN, mprice, gTol, subTol )
%
% Computes reduced gradients for variables listed in set N
% and returns up to mprice most significant ones dT
% (bigger than optTol) and their positions kT in N.
% Their positions in A are given by T = N(kT).
%
% 05 Sep 2006: First version of lsl1price.
%              Full pricing.  Finds just one variable T = N(kT).
% 06 Sep 2006: First version of multiple price.
% 09 Nov 2007: Renamed from lsl1price2 to BPprice. (Used in BPprimal.)
% 09 Aug 2008: delta is a new input parameter.
%              Objective is now ||x||_1 + 1/2 delta x'x + 1/2 lambda y'y.
%              Reduced gradients for x are     (sign(x) + delta x) - A'y.
% 18 Aug 2008: delta is irrelevant if we assume xN = 0.

% AUTHORS
% Michael P. Friedlander               Michael A. Saunders
% University of British Columbia       Stanford University
% mpf@cs.ubc.ca                        saunders@stanford.edu
%
% $Id: BPprice.m 385 2008-09-02 22:22:41Z mpf $
%----------------------------------------------------------------------

  dT     =  [];
  kT     =  [];
  rgmax  =  0;
  Nty    =  Aty(iN);
  lenN   =  length(Nty);

  nneg   = x(iN) >= 0;
  npos   = x(iN) <= 0;
  
  % Reduced gradients:
  % xN increasing: we want negative ones
  g1     = - Nty + nneg - ~nneg; % -Nty + 1 if x >= 0; -Nty - 1 if x < 0

  % xN decreasing: we want positive ones
  g2     = - Nty - npos + ~npos; % -Nty - 1 if x <= 0; -Nty + 1 if x > 0

  g12    =  [g1; -g2];
  [gs,p] =  sort(g12);       % g12 sorted in ascending order

  lenT   =  min(mprice,lenN);% Keep up to mprice best ones
  if lenT == 0, return, end

  rgmax  =  abs(gs(1));      % Always >= 0
  if gs(1) >= -gTol
     rgmax  = 0;
     dT     = [];
     kT     = [];
  else
     gs     =  gs(1:lenT);
     sTol   =  max( subTol*rgmax, gTol );  % Keep only significant ones
     lenT   =  sum( gs < -sTol );
     if lenT > 0
        dT     =  gs(1:lenT);
        kT     =  p(1:lenT);
        I2     =  find(kT > lenN);      % fix entries from g2
        dT(I2) = -dT(I2);
        kT(I2) =  kT(I2) - lenN;
     end
  end
end % function BPprice
