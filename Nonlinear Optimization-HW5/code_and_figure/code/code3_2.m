clear;
clc;
load('HW5_data.mat');
mu = 1e0;
A = Q + mu*c*c';

n = length(c);
x_0 = zeros(n,1);

M = 10000;
x = zeros(n,M);
x(:,1) = x_0;
error = 1e-4 ;
number_itegrate = 0;
grad = A*x(:,number_itegrate+1)-b;

while norm(grad) >= error
    number_itegrate = number_itegrate + 1;

    alpha = 1;
    eta  = 1.1;
    varepsion = 0.2;
    d = -grad;
    x_1 = x(:,number_itegrate);
    x_2 = x_1 + alpha*d;
  while (1/2*x_2'*A*x_2 - b'*x_2) - ...
            (1/2*x_1'*A*x_1 - b'*x_1) > -varepsion*alpha*(d'*d)  
      alpha = alpha/eta;
      x_2 = x_1 + alpha*d;
  end
  
    x(:,number_itegrate+1) = x(:,number_itegrate)  - alpha*grad;
    grad = A*x(:,number_itegrate+1)-b;

end


n = length(c);
x_0 = zeros(n,1);
x_star =  gaussseidel(Q+mu*c*c',b,x_0,1e-9,3e5);




X = 0:1:number_itegrate;
Y1 = zeros(number_itegrate+1,1);
Y2 = zeros(number_itegrate+1,1);
Y3 = zeros(number_itegrate+1,1);
for i = 1:number_itegrate+1
  Y1(i) = norm(A*x(:,i)-b);
  Y2(i) = 1/2*x(:,i)'*Q*x(:,i) - b'*x(:,i) + mu/2*(c'*x(:,i))^2;
  Y3(i) = norm(x(:,i)-x_star);
 end

figure(1)
semilogy(X,Y1,'LineWidth',2);
xlabel('k');
ylabel('||\nabla f(x_k)||');
saveas(gca,'321_1.eps');


figure(2)
plot(X,Y2,'LineWidth',2);
xlabel('k');
ylabel('f(x_k)');
saveas(gca,'321_2.eps');

figure(3)
semilogy(X,Y3,'LineWidth',2);
xlabel('k');
ylabel('||x_k-x*||');
saveas(gca,'321_3.eps');

