function R = QRaddrow(R,a)
%QRADDROW   Add a row and update a Q-less QR factorization.
%
% R1 = QRaddrow(R,a) returns the triangular part of a QR factorization
% of [A; a], where A = QR for some Q.  The argument 'a' should be a
% row vector.

% 05 Apr 2008: First version of QRaddrow.
%              Michael Friedlander (mpf@cs.ubc.ca) and
%              Michael Saunders (saunders@stanford.edu)

  n = size(R,1);

  for k=1:n
      G = planerot( [R(k,k); a(k)] );
      
      B = G * [ R(k,k:n)
                a(  k:n) ];
      
      R(k,k:n) = B(1,:);
      a(  k:n) = B(2,:);
  end
