
% -*- texinfo -*-
% @deftypefn  {} {[@var{mu}, @var{sigma}, @var{f}] =} util_lanczos (@var{A}, @
% @var{At}, @var{tau}, @var{xi}, @var{y}, @var{y0})
% Compute the solution to a linear system and the diagonal of the
% inverse coefficient matrix using Lanczos iteration.
%
% The outputs @var{mu} and @var{sigma} are the solution and
% diagonal, and @var{f} holds an estimate of the determinant
% of the kernel matrix.
% @end deftypefn
%
function [mu, sigma, f] = util_lanczos (A, At, tau, xi, y, y0)
  % get the problem sizes.
  n = length(xi);
  m = length(y);

  % define the kernel matrix operator.
  K = @(z) z ./ tau + A(At(z) ./ xi);

  % initialize the tridiagonal matrix elements.
  omega = [];
  gamma = [];

  % initialize the lanczos vector.
  q = y ./ norm(y);
  Q = q;

  % initialize the transformed lanczos vector.
  v = At(q);
  V = v;

  % compute the first conjugate vector.
  p = K(q);
  omega = real(q' * p);
  p = p - omega .* q;

  % initialize the recurrence variables.
  f2 = 0;
  f1 = 1;
  f = omega;

  % run lanczos iterations.
  for j = 2 : m
    % fully reorthogonalize.
    h = Q(:, 1 : j-1)' * p;
    p = p - Q(:, 1 : j-1) * h;

    % compute the next lanczos vector.
    gamma = [gamma; norm(p)];
    q = p ./ gamma(end);
    Q = [Q, q];

    % compute the next transformed lanczos vector.
    v = At(q);
    V = [V, v];

    % compute the next conjugate vector.
    p = K(q);
    omega = [omega; real(q' * p)];
    p = p - omega(end) .* q - gamma(end) .* Q(:, j-1);

    % update the recurrence variables.
    f2 = f1;
    f1 = f;
    f = omega(end) * f1 - gamma(end)^2 * f2;
  end

  % construct the tridiagonal matrix.
  T = sparse(diag(omega) + diag(gamma, 1) + diag(gamma, -1));

  % compute the marginal means.
  Qty = sparse(norm(y) .* eye(m, 1));
  t = Q * (T \ Qty);
  mu = At(t) ./ xi;

  % compute the marginal variances.
  S = sum(real(V .* (T \ V').'), 2);
  sigma = 1./xi - S ./ xi.^2;
end

