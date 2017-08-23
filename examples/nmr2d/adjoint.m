
% define the adjoint operator.
function x = adjoint (y, schd, n)
  % create an infilled vector of y-values.
  f = zeros(n(1), n(2), 4);
  f(schd) = y;

  % compute the fourier transform.
  x = mxfft2(f, n);
end

