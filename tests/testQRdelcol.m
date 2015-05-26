function testQRdelcol(n)

%        testQRdelcol(n)
% generates a random nxn A, scales it to be ill-conditioned,
% computes R = qr(A) and updates it as each column of R is deleted.

% 18 Jun 2007: First version of testQRdelcol.m.

rand('state',0)
A  = diag(1./(1:n))*rand(n,n);
R  = triu(qr(A));
nR = n;

for j=2:n
    R      = QRdelcol(R,2);  % Delete the 2nd column each time
    A(:,2) = [];             % Delete column of A
    nR     = nR - 1;
    R2     = triu(qr(A));
    R2     = R2(1:nR,1:nR);
    E      = abs(R2) - abs(R);
    E2     = E.*E;
    emax   = sqrt(max(sum(E2)))
    dmin   = min(abs(diag(R)))
end

return
