function testQRaddcol(n)

%        testQRaddcol(n)
% generates a random nxn A, scales it to be ill-conditioned,
% and builds up R column-wise such that R'R = A'A.

% 15 Jun 2007: First version of testQRaddcol.m.
%              Note that each column of R above the diagonal
%              could differ in sign from QR's R, even if the
%              diagonal elements agree.
% 18 Jun 2007: R now has the exact size each time.

rand('state',0)
A = diag(1./(1:n))*rand(n,n);

% Test 1: beta = 0
fprintf('Test 1: beta = 0\n')

R    = [];
S    = zeros(n,0);

for j=1:n
    a   = A(:,j);
    R   = QRaddcol( S, R, a );
    S   = [S a];
    SS  = S'*S;
    err = norm( SS - R'*R,'fro');
    fprintf('%2i %15.6e\n',j,err)
end

% Test 2: beta > 0
beta = 1e-1;
fprintf('Test 2: beta = %g\n',beta)
R    = [];
S    = zeros(n,0);

for j=1:n
    a   = A(:,j);
    R   = QRaddcol( S, R, a, beta );
    S   = [S a];
    SS  = S'*S;
    err = norm( SS + beta^2*eye(j) - R'*R,'fro');
    fprintf('%2i %15.6e\n',j,err)
end
