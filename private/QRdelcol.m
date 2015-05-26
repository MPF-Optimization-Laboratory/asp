function R = QRdelcol(R,k)
%QRDELCOL   Delete the k-th column and update a Q-less QR factorization.
%
% R = QRdelcol(R,k) deletes the k-th column of upper-triangular R and
% restores it to upper-triangular form.  On input, R is an n x n upper
% triangular matrix.  On output, R is an n-1 x n-1 upper triangle.

% 18 Jun 2007: First version of QRdelcol.m.
%              Michael Friedlander (mpf@cs.ubc.ca) and
%              Michael Saunders (saunders@stanford.edu)
%              To allow for R being sparse,
%              we eliminate the k-th row of the usual
%              Hessenberg matrix rather than its subdiagonals,
%              as in Reid's Bartel-Golub LU update and also
%              the Forrest-Tomlin update.
%              (But Matlab will still be pretty inefficient.)
% 18 Jun 2007: R is now the exact size on entry and exit.

  R(:,k) = [];      % Delete the k-th column
  [~,n]  = size(R); % Should have m=n+1

  for j=k:n         % Forward sweep to reduce k-th row to zeros
    I = [j+1 k];
    [G,y]    = planerot(R(I,j));
    R(j+1,j) = y(1);
    if j<n && G(2,1)~=0
      R(I,j+1:n) = G*R(I,j+1:n);
    end
  end

  R(k,:) = [];      % Delete the k-th row
  return
