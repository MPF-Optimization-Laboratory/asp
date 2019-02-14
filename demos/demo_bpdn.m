% Setup a simple test problem.
randn('state',0); rand('state',0);
m    = 600; n = 2560; k = 20;      % No. of rows, columns, and nonzeros
p    = randperm(n);   p = p(1:k);  % Position of nonzeros in x
x    = zeros(n,1);                 % Generate sparse solution
x(p) = randn(k,1);
A    = randn(m,n);                 % Gaussian m-by-n ensemble
b    = A*x;                        % Compute the RHS vector

% Choose a sequence of 10 linearly-spaced values of lambda.
lambdamax = norm(A'*b,inf);
lambdas   = linspace(lambdamax,0,10);

opts          = as_setparms;
opts.loglevel = 1;
inform        = [];  % IMPORTANT: must initialize in this way.

% Re-solve the BPDN problem: minimize .5||Ax-b||^2 + lambda||x||_1
for lambda = lambdas
    [x,inform] = as_bpdn(A,b,lambda,opts,inform);
    
    % x is now optimal for the current lambda.
    % Here is where you can examine the solution and perhaps exit the loop.
    stem(x);
    xlim([0,n]);
    ylim([-2.5,+2.5]);
    pause(0.05);
end


