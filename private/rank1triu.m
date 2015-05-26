function R = rank1triu(R, u, v)
% rank1triu  triangularize Q*(R + u v'), where R is triangular; discard Q
%
% R = rank1triu(R,u,v)

% 12 Aug 07: First version of rank1triu.
%% $Id$

[m,n] = size(R);  % R should be square
if n==0
    return
end

G = zeros(2,2); % Preallocate Givens matrix.

% Backward sweep to reduce u(1:n-1) to zeros.
for j = n-1:-1:1
    if u(j) == 0, cycle; end
    I      = [n j];
    r      = norm(u(I));
    G(1,1) =  u(n);
    G(1,2) =  u(j);
    G(2,1) = -u(j);
    G(2,2) =  u(n);
    G      = G / r;
    u(n)   = r;
    R(I,j:n) = G*R(I,j:n);
end
R(n,:) = R(n,:) + u(n)*v';

% R now has a spike in its last row.
% Forward sweep to restore triangular shape.
for j = 1:n-1
    if R(n,j) == 0, cycle; end
    I          = [j n];
    r          = norm(R(I,j));
    G(1,1)     =  R(j,j);
    G(1,2)     =  R(n,j);
    G(2,1)     = -R(n,j);
    G(2,2)     =  R(j,j);
    G          = G / r;
    R(j,j)     = r;
    R(n,j)     = 0;
    R(I,j+1:n) = G*R(I,j+1:n);
end

end % function rank1triu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function test
%% To test the code:
n = 100;
kappa = 1e8;
randn('state',0); rand('state',0);
A = gallery('randsvd',n,kappa,3);
u = randn(n,1);
v = randn(n,1);
[Q,R] = qr(A);
Rb1 = rank1triu(R,u,v);
Rb2 = R + u*v';
err = norm(Rb1'*Rb1 - Rb2'*Rb2);
fprintf('%15.6e\n', err);
end
