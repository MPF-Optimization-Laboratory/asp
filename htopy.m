m = 600; n = 2560; k = 20;    % No. of rows, columns, and nonzeros
p = randperm(n); p = p(1:k);  % Position of nonzeros in x
x = zeros(n,1);               % Generate sparse solution
x(p) = randn(k,1);
A = randn(m,n);               % Gaussian m-by-n ensemble
b = A*x;                      % Compute the RHS vector
x = as_ktopy(A,b,10);         % Find a 10-sparse estimate of the soln
