
% set the problem sizes.
m = 200;
n = 1000;
k = 5;

% set the prior parameters.
mu0 = 200;
lambda0 = 0.02;
alpha0 = 1e-9;
beta0 = 1e-9;

% set the iteration count.
iters = 10000;

% set the smoothing kernel width.
delta = 50;

% build a measurement matrix.
A = exp(-([1 : n/m : n]' - [1 : n]).^2 ./ delta^2);
A = diag(1 ./ sqrt(sumsq(A, 2))) * A;
At = [];

% build the ground truth vector.
id0 = [100 : n - 100]';
idx = id0(randperm(length(id0))(1 : k));
x0 = zeros(n, 1);
x0(idx) = 1;

% sample a noise vector.
z = normrnd(0, 0.01, m, 1);

% construct the measurement.
y = A * x0 + z;

% save the current instance.
save('-binary', '-z', 'instance.dat.gz');

