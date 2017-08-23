
% define the forward operator.
function y = forward (x, schd, n)
  % inverse fourier transform and subsample.
  f = mxifft2(x, n);
  y = f(schd);
end

