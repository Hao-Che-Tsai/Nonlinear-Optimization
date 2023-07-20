clc, clear
Q = [2 1 0 10;
   1 4 3 0.5;
0 3 -5 6;
   10 0.5 6 -7];
b = [-1; 0; -2; 3];
vareplison = 1e-9;
beta = 2;
c = 1;
k = 1;
F = zeros(80,1);
x = zeros(4,80);
p = zeros(80,1);
cP = zeros(80,1);
f = @(x) 1/2*x'*Q*x + b'*x;
P = @(x) 1/2*((sum((max(0,-x)).^2)) + (sum(x.^2)-1)^2);
func = @(x) 1/2*x'*Q*x + b'*x + c/2*((sum((max(0,-x)).^2)) + (sum(x.^2)-1)^2);
x(:,k) = fminsearch(func,zeros(4,1));
while P(x(:,k)) > vareplison
   k = k+1;
c = c*beta;
func = @(x) 1/2*x'*Q*x + b'*x + c/2*((sum((max(0,-x)).^2)) + (sum(x.^2)-1)^2);
   x(:,k) = fminsearch(func,x(:,k-1));
end
for i = 1:k
   F(i) = f(x(:,i));
   p(i) = P(x(:,i));
   cP(i) = 2^(i-1)*p(i);
end
xfinal = x(:,k)
minf = f(xfinal)
A = zeros(3,4);
A(1,2) = -1;
A(2,4) = -1;
A(3,:) = xfinal'*2;
nabla = Q*xfinal + b;
lambda = (A*A')^(-1)*A*(-nabla)
r = null(A);
L = Q+2*lambda(3)*eye(4)
r'*L*r
t = 1:k;
figure(1)
plot(t,F(1:k,1)');
xlabel('k');
ylabel('f(x_k)');
figure(2)
semilogy(t,cP(1:k,1)');
xlabel('k');
ylabel('c_kP(x_k)');
figure(3)
semilogy(t,p(1:k,1)');
xlabel('k');
ylabel('P(x_k)');