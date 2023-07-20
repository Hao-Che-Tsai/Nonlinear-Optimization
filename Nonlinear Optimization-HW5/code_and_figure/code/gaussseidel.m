function [x,k] = gaussseidel(A,b,x_0,ep,N)

n = length(b);
k = 0;
if nargin < 5
    N = 1e3;
end
if nargin < 4
    ep = 1e-5;
end
if nargin < 3
    x_0 = zeros(n,1);
 
end

x = x_0;
x_0 = x + 2*ep;

U1 = tril(A);
U2 = U1^(-1);
G = -U2*(A-U1); f = U2*b;

 while norm(x_0-x,inf) > ep && k < N
     k = k + 1;
     x_0 = x;
     x = G*x_0 + f;
 end

