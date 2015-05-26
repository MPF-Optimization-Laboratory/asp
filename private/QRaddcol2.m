function R = QRaddcol(A,R,a)

%        R = QRaddcol(A,R,a);
% adds a column to the QR factorization of A without using Q.
% On entry, R is a DENSE n x n upper triangular matrix such that A = Q*R.
% A should have fewer columns than rows.
% A and a may be sparse or dense.
% On exit, R is the n+1 x n+1 factor corresponing to [A a] = Qnew*Rnew.
%
%   Anew = [A  a],   Rnew = [R   u  ],   Rnew'Rnew = Anew'Anew.
%                           [  gamma]
% The new column "a" is assumed to be nonzero.
% If A has no columns yet, input A = [], R = [];

% 15 Jun 2007: First version of QRaddcol.m.
%              Where necessary, Ake Bjorck's CSNE
%              (Corrected Semi-Normal Equations) method
%              is used to improve the accuracy of u and gamma.
%              See p143 of Ake's Least Squares book.
% 18 Jun 2007: R is now the exact size on entry and exit.
% 19 Oct 2007: Sparse A, a makes c and u sparse.
%              Force them to be dense.
%              For dense R we probably should use linsolve,
%              which requires c and u to be dense anyway.
% 04 Aug 2008: This is QRaddcol2.m.
%              It updates u using du, rather than u = R*z as in Ake's book.
%              We guess that it might be slightly more accurate,
%              but it's hard to tell.  No R*z makes it a little cheaper.
%
% $Id: QRaddcol.m 351 2008-08-05 05:45:24Z saunders $

  [m,n]  = size(R);
% if m~=n, error('R must be square')
  anorm  = norm(a);
  if n==0
    R = anorm;
    return
  end

  c      = A'*a;
  u      = full(R'\c);    % In case sparse A,a,c made u sparse.
  anorm2 = anorm^2;
  unorm2 = norm(u)^2;
  d2     = anorm2 - unorm2;

  if d2 > anorm2 %DISABLE 0.01*anorm2     % Cheap case: gamma is not too small
    gamma = sqrt(d2);
  else
    z     = R\u;          % First approximate solution to min ||Az - a||
    r     = a - A*z;
    c     = full(A'*r);
    du    = R'\c;
    dz    = R\du;
    u     = u + du;       % Refinement of both u and z
    z     = z + dz;
    r     = a - A*z;
    gamma = norm(r);      % Safe computation (we know gamma >= 0).
  end

  R = [    R         u
       zeros(1,n)  gamma];
  return
