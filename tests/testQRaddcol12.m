function testQRaddcol12(n,kappa)

%        testQRaddcol12(n,kappa)
%        testQRaddcol12(50,1e12)
%
% Compares two versions of Ake Bjorck's CSNE approach
% to adding a column to R in a QR factorization without Q.
% Uses gallery randsvd matrix A of order n with cond(A) = kappa.
% QRaddcol0 follows Ake's version in using u = R*z.
% QRaddcol  refines both u and z.
% Sometimes there is a significant difference in the accuracy of R,
% but not always favoring one method.
% The modified QRaddcol is slightly cheaper because it doesn't need R*z.

% 04 Aug 2008: First version of testQRaddcol.  MAS and MPF.

randn('state',0); rand('state',0);
A  = gallery('randsvd',n,kappa,3);  % 3 means geometrically distributed sv's.
S  = zeros(n,0);
R1 = [];
R2 = [];

for j=1:n
  a    = A(:,j);
  R1   = QRaddcol1(S,R1,a);
  R2   = QRaddcol2(S,R2,a);
  S    = [S a];
  SS   = S'*S;
  err1 = norm(SS - R1'*R1,'fro');
  err2 = norm(SS - R2'*R2,'fro');
  fprintf('%15.6e %15.6e\n',err1, err2)
end

end % function testQRaddcol
