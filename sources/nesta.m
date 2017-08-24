
% nesterov's algorithm for quadratically constrained l1-minimization.
%
function [x, elbo, eta] = ...
nesta (b, A, At, mu0, lambda0, alpha0, beta0, iters)
  % adjust the iteration scheme to fit nesta.
  if (isscalar(iters))
    iters = [round(iters / 50), 50];
  end

  % create a function handle for projection.
  AtA = @(x) At(A(x));

  % define the gradient of the l1 term.
  df = @(x, mu) (x ./ mu) .* (abs(x) < mu) + sign(x) .* (abs(x) >= mu);

  % initialize the transformed data vector.
  h = At(b);

  % get the problem sizes.
  m = length(b);
  n = length(h);

  % initialize the iterates.
  x = zeros(n, 1);
  y = zeros(n, 1);
  z = zeros(n, 1);

  % set up the vector of thresholds, and compute the max threshold.
  muv = linspace(0.9, 0.01, iters(1));
  mu0 = max(abs(h));

  % compute the constant for the quadratic constraint.
  vareps = 0.1 * sqrt(m / (mu0 / lambda0));

  % iterate, outer loop.
  for it = 1 : iters(1)
    % set the current threshold value and lipschitz constant.
    mu = mu0 * muv(it);
    L = 1 / mu;

    % initialize the gradient accumulator.
    xc = x;
    adf = zeros(size(x));

    % iterate, inner loop.
    for jt = 1 : iters(2)
      % compute the gradient.
      dfdx = df(x, mu);

      % compute the scale factors for the accumulator and iterates.
      alpha = (jt + 1) / 2;
      tau = 2 / (jt + 3);

      % update the accumulator.
      adf += alpha .* dfdx;

      % update the y-iterate.
      ly = max(0, (L/vareps) * norm(b - A(x - dfdx ./ L)) - L);
      y = x + (ly / L) .* h - dfdx ./ L;
      y -= (ly / (ly + L)) .* AtA(y);

      % update the z-iterate.
      lz = max(0, (L/vareps) * norm(b - A(xc - adf ./ L)) - L);
      z = xc + (lz / L) .* h - adf ./ L;
      z -= (lz / (lz + L)) .* AtA(z);

      % update the x-iterate.
      x = tau .* z + (1 - tau) .* y;
    end
  end

  % return empty values for the other outputs.
  elbo = [];
  eta = [];
end
