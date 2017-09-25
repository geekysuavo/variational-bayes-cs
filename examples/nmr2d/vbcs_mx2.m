
% modification of vbcs() for two-dimensional multicomplex data,
% under the assumption of a partial fourier measurement.
%
function [x, elbo, eta] = ...
vbcs_mx2 (y, A, At, mu0, lambda0, alpha0, beta0, iters)
  % check for the minimum number of arguments.
  if (nargin < 3)
    error('at least three arguments required');
  end

  % check the required measurement vector.
  if (isempty(y) || !ismatrix(y) || columns(y) != 4)
    error('measurement must be a 4-column matrix');
  end

  % check the required function handles.
  if (isempty(A) || !is_function_handle(A) || ...
      isempty(At) || !is_function_handle(At))
    error('invalid measurement: expected function handles');
  end

  % create a function handle for projection.
  AtA = @(x) At(A(x));

  % check for a prior mu parameter.
  if (nargin < 4 || isempty(mu0))
    % none specified. use a default value.
    mu0 = 1e-3;
  end

  % check for a prior lambda parameter.
  if (nargin < 5 || isempty(lambda0))
    % none specified. use a default value.
    lambda0 = 1e-3;
  end

  % check for a prior alpha parameter.
  if (nargin < 6 || isempty(alpha0))
    % none specified. use a default value.
    alpha0 = 1e-3;
  end

  % check for a prior beta parameter.
  if (nargin < 7 || isempty(beta0))
    % none specified. use a default value.
    beta0 = 1e-3;
  end

  % check for an iteration count argument.
  if (nargin < 8 || isempty(iters))
    % none specified. use a default value.
    iters = 100;
  end

  % initialize the transformed data vector.
  h = At(y);

  % get the problem sizes.
  m = rows(y);
  n = prod(size(h)(1:2));

  % initialize the weight means.
  x = zeros(size(h));

  % initialize the precisions.
  xi = repmat(alpha0 ./ beta0, size(h,1), size(h,2));

  % initialize the noise.
  tau = mu0 / lambda0;

  % initialize the projected weight means.
  u = AtA(x);

  % compute the projector diagonal, valid only for partial fourier.
  g = m / n;

  % compute the constant updated parameters.
  alpha = alpha0 + (1/2);
  mu = mu0 + m/2;

  % compute the constant portion of the lower bound.
  L0 = n * alpha0 * log(beta0) + mu0 * log(lambda0) ...
     + n * (lgamma(alpha) - lgamma(alpha0)) ...
     + (lgamma(mu) - lgamma(mu0)) ...
     + n * alpha + mu - (m/2) * log(2*pi) + n/2;

  % initialize the lower bound vector.
  elbo = repmat(L0, iters, 1);

  % compute the inner product of the measurement.
  yy = mxdot(y, y)(1);

  % iterate.
  for it = 1 : iters
    % update the variances.
    v = 1 ./ (xi + tau .* g);

    % compute the parallel means and their projection.
    xp = tau .* bsxfun(@times, v, h - u + g .* x);
    up = AtA(xp);

    % compute the terms of the step size.
    T1 = tau .* mxdot2(h, x - xp);
    T2 = -mxdot2(x,  tau .* u  + bsxfun(@times, xi, x));
    T3 = -mxdot2(xp, tau .* up + bsxfun(@times, xi, xp));
    T4 =  mxdot2(x,  tau .* up + bsxfun(@times, xi, xp));

    % compute the bounded step size for the mean update.
    gamma = (T1 + T2 + T4)(1) / (T2 + T3 + 2*T4)(1);
    gamma = min(max(gamma, 1e-3), 1);

    % update the means and their projection.
    x = (1 - gamma) .* x + gamma .* xp;
    u = (1 - gamma) .* u + gamma .* up;

    % compute some intermediate values.
    trGxx = mxdot2(x, u)(1) + sum(vec(g .* v));
    hx = mxdot2(h, x)(1);
    x2 = sumsq(x, 3) + v;

    % update the precisions.
    beta = beta0 + 0.5 .* x2;
    xi = alpha ./ beta;

    % update the noise.
    lambda = lambda0 + 0.5 .* yy - hx + 0.5 .* trGxx;
    lambda = real(lambda);
    tau = mu / lambda;

    % compute the new value of the lower bound.
    L = -(tau/2) * (yy - 2 * hx + trGxx) ...
      - sum(vec(((alpha/2) * x2 + alpha * beta0) ./ beta)) ...
      - mu * log(lambda) - mu * lambda0 / lambda ...
      - alpha * sum(vec(log(beta))) ...
      + 0.5 * sum(vec(log(v)));

    % store the bound.
    elbo(it) += L;
  end

  % check if the complete set of variational parameters was requested.
  if (nargout >= 3)
    % store the parameters over x.
    eta.xhat = x;
    eta.v = v;

    % store the parameters over xi.
    eta.alpha = alpha;
    eta.beta = beta;

    % store the parameters over tau.
    eta.mu = mu;
    eta.lambda = lambda;
  end
end

