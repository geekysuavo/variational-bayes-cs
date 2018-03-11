
% set the problem sizes.
m = 200;
n = 1000;
k = 20;

% set the prior parameters.
nu0 = 200;
lambda0 = 0.005;
alpha0 = 1e-3;
beta0 = 1e-3;

% set the iteration count.
iters = 500;

% build a measurement matrix.
A = normrnd(0, 1, m, n);
A = diag(1 ./ sqrt(sumsq(A, 2))) * A;
At = [];

% sample properties for a ground truth vector.
sgn = 2 * (randi(2, k, 1) - 1.5);
amp = unifrnd(0.1, 1.1, k, 1);
idx = randperm(n)(1 : k)';

% build the ground truth vector.
x0 = zeros(n, 1);
x0(idx) = sgn .* amp;

% sample a noise vector.
z = normrnd(0, 0.005, m, 1);

% construct the measurement.
y = A * x0 + z;

% save the current instance.
save('-binary', '-z', 'instance.dat.gz');

