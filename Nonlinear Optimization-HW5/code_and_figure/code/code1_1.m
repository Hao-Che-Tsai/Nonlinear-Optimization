x_0 = zeros(2,1);
M = 1000;
x = zeros(2,M);
x(:,1) = x_0;

epsion = 0;
eta  = 1.1;
varepsion = 0.5;
delta = 1;
error = 1e-10;
number_itegrate = 0;


while norm(grad(x(:,number_itegrate+1))) > error
  number_itegrate = number_itegrate + 1;
  
  
  x_1 = x(:,number_itegrate);
  
  hess = hessian(x_1);
  epsion = max(delta - min(eig(hess)),0);
  g = grad(x_1);
  d = -(epsion*eye(2)+hess)^(-1)*grad(x_1);
  alpha = 1;
  x_2 = x_1 + alpha*d;

  while f(x_2) - f(x_1) > -varepsion*alpha*(g'*d) 
      alpha = alpha/eta;
      x_2 = x_1 + alpha*d;
    end
  
  
   x(:,number_itegrate+1) = x(:,number_itegrate)  + alpha*d;
end



x_star = x(:,number_itegrate+1);

norm(grad(x_star))
hessian(x_star)

X = 0:1:number_itegrate;
Y1 = zeros(number_itegrate+1,1);
Y2 = zeros(number_itegrate+1,1);
for i = 1:number_itegrate+1
  Y1(i) = norm(grad(x(:,i)));
  Y2(i) = f(x(:,i));
 end

figure(1)
semilogy(X,Y1,'LineWidth',2);
xlabel('k');
ylabel('||\nabla f(x_k)||');
saveas(gca,'1_1.eps');


figure(2)
plot(X,Y2,'LineWidth',2);
xlabel('k');
ylabel('f(x_k)');
saveas(gca,'1_2.eps');

