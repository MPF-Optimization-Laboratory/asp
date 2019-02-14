function tests = testASP
tests = functiontests(localfunctions);
end

function test_asp_bpdn(testCase)
m = 600; n = 2560; k = 20;    % No. of rows, columns, and nonzeros
p = randperm(n); p = p(1:k);  % Position of nonzeros in x
x = zeros(n,1);               % Generate sparse solution
x(p) = randn(k,1);
A = randn(m,n);               % Gaussian m-by-n ensemble
b = A*x;                      % Compute the RHS vector
opts.loglevel = 0;
opts = as_setparms(opts);
[xs, inform] = as_bpdn(A,b,0,opts);  % Solve the basis pursuit problem
y = inform.y;

pFeas = norm(A*xs - b,Inf)/max(1,norm(xs,Inf));
dFeas = max(0,norm(A'*y,Inf) - 1);
dComp = abs(norm(x,1) - b'*y);

verifyLessThanOrEqual(testCase, pFeas, 1e-6);
verifyLessThanOrEqual(testCase, dFeas, 1e-6);
verifyLessThanOrEqual(testCase, dComp, 1e-6);


end
