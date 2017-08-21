
% load the input data.
load('fid.dat.gz');
x0 = x0 ./ max(real(x0));

% =========================================================

% define the forward operator.
function y = forward (x, schd, n)
  % inverse fourier transform and subsample.
  f = ifft(x) .* sqrt(n);
  y = f(schd,:);
end

% define the adjoint operator.
function x = adjoint (y, schd, n)
  % create an infilled vector of y-values.
  f = zeros(n, columns(y));
  f(schd,:) = y;

  % compute the fourier transform.
  x = fft(f) ./ sqrt(n);
end

% =========================================================

% set the problem sizes.
m = length(sched);
n = 2 * length(x0);

% set the prior parameters.
mu0 = 50;
lambda0 = 0.005;
alpha0 = 1e-3;
beta0 = 1e-3;

% set the iteration count.
iters = 500;

% define the final operators.
A = @(x) forward(x, sched, n);
At = @(y) adjoint(y, sched, n);

% construct the subsampled measurement.
y = x0(sched);

% save the current instance.
save('-binary', '-z', 'instance.dat.gz');

