function R = QRaddcol(A,R,a,beta)
%QRADDCOL   Add a column to a QR factorization without using Q.
%
% R = QRaddcol(A,R,V)  adds the m-vector V to the QR factorization of the
% m-by-n matrix A without using Q. On entry, R is a dense n-by-n upper
% triangular matrix such that R'*R = A'*A.
%
% R = QRaddcol(A,R,V,beta)  is similar to the above, except that the
% routine updates the QR factorization of [A; beta I], and R'*R= (A'*A +
% beta^2*I) = R'*R.
%
% A should have fewer columns than rows.  A and a may be sparse or dense.
% On exit, R is the (n+1)-by-(n+1) factor corresponding to
%
%   Anew = [A        V ],    Rnew = [R   u  ],   Rnew'Rnew = Anew'Anew.
%          [beta*I     ]            [  gamma]
%          [       beta]
%
% The new column V is assumed to be nonzero.
% If A has no columns yet, input A = [], R = [];

% 15 Jun 2007: First version of QRaddcol.m (without beta).
%              Michael Friedlander (mpf@cs.ubc.ca) and
%              Michael Saunders (saunders@stanford.edu)
%              Where necessary, Ake Bjorck's CSNE
%              (Corrected Semi-Normal Equations) method
%              is used to improve the accuracy of u and gamma.
%              See p143 of Ake's Least Squares book.
% 18 Jun 2007: R is now the exact size on entry and exit.
% 19 Oct 2007: Sparse A, a makes c and u sparse.
%              Force them to be dense.
%              For dense R we probably should use linsolve,
%              which requires c and u to be dense anyway.
% 04 Aug 2008: QRaddcol2 updates u using du, rather than
%              u = R*z as in Ake's book.
%              We guess that it might be slightly more accurate,
%              but it's hard to tell.  No R*z makes it a little cheaper.
% 03 Sep 2008: Generalize A to be [A; beta*I] for some scalar beta.
%              Update u using du, but keep Ake's version in comments.

[m,n]  = size(R);
% if m~=n, error('R must be square')
if nargin < 4 || isempty(beta)
   beta = 0;
end
anorm  = norm(a);
anorm2 = anorm^2;
beta2  = beta^2;
if beta~=0
   anorm2 = anorm2 + beta2;
   anorm  = sqrt(anorm2);
end

if n==0
   R = anorm;
   return
end

c      = A'*a;          % = [A' beta*I 0]*[a; 0; beta]
u      = full(R'\c);    % In case sparse A,a,c made u sparse.
unorm2 = norm(u)^2;
d2     = anorm2 - unorm2;

if d2 > anorm2 %DISABLE 0.01*anorm2     % Cheap case: gamma is not too small
   gamma = sqrt(d2);
else
   z     = R\u;          % First approximate solution to min ||Az - a||
   r     = a - A*z;
   c     = full(A'*r);
   if beta~=0,
      c = c - beta2*z;
   end
   du    = R'\c;
   dz    = R\du;
   z     = z + dz;       % Refine z
   % u     = R*z;          % Original:     Ake's version.
   u     = u + du;       % Modification: Refine u
   r     = a - A*z;
   gamma = norm(r);      % Safe computation (we know gamma >= 0).
   if beta~=0
      gamma = sqrt(gamma^2 + beta2*norm(z)^2 + beta2);
   end
end

R = [    R         u
     zeros(1,n)  gamma];
