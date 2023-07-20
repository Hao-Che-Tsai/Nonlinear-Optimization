function y = hessian(x)
  x1 = x(1,1);
  x2 = x(2,1);
  y(1,1) = 2;
  y(1,2) = -5;
  y(2,1) = -5;
  y(2,2) = 12*x2^2;