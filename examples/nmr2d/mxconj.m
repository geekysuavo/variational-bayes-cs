
function y = mxconj (x)
  y = x;
  y(:,1) *= -1;
  y(:,2) *= -1;
end

