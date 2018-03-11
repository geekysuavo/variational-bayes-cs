
% -*- texinfo -*-
% @deftypefn {} {@var{phi0} =} vrvm_const (@var{m}, @var{n}, @
% @var{nu0}, @var{lambda0}, @var{alpha0}, @var{beta0})
% Compute the constant offset term of the Variational Relevance
% Vector Machine objective function, based only on prior values
% and the problem size.
% @end deftypefn
%
function phi0 = vrvm_const (m, n, nu0, lambda0, alpha0, beta0)
  % compute Cp, the constant offset of the Kullback-Liebler corrected
  % variational lower bound.
  Cp = nu0 * log(lambda0) + n * alpha0 * log(beta0) ...
     + (lgamma(nu0 + m/2) - lgamma(nu0)) ...
     + n * (lgamma(alpha0 + 1/2) - lgamma(alpha0)) ...
     - ((m + n) / 2) * log(2 * pi);

  % compute phi0, the constant offset of the objective function,
  % which includes the constant term of the normal entropy.
  phi0 = -(n/2) * log(2 * pi * e) - Cp;
end

