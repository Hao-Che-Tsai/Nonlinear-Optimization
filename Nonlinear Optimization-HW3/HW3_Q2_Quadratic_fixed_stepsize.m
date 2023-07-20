% Implementation of Gradient Descent with a fixed step size for Quadratic
% optimization problem
% minimize f(x) = 1/2 x^T Q x + a^{T} x


Q = [1 0.5 0;
    0.5 1 0.25;
    0 0.25 1;]; % Q is positive definite.

a = [1; -1; 2];

x_opt = Q\(-a);

lambda_max = max(eig(Q)); % max eigenvalue of lambda

alpha = 1/lambda_max % The fixed step size alpha should be smaller that 2/lambda_max to ensure global convergence

%% Gradient Descent

itr = 1;
x = [0; 0; 0]; % initial point

epsilon = 10^(-6);

while (norm(Q*x(:,itr) + a,2) > epsilon) % Stopping criterion: Euclidean norm of gradient <= epsilon 


    x_new = x(:,itr) - alpha*(Q*x(:,itr) + a); % Gradient Descent step
    x = [x x_new];
    itr = itr + 1;
    itr
end

% grad_x will save the Euclidean norm of gradient at each iteration.
grad_x = [];

for k=1:itr
   
    grad_x = [grad_x; norm(Q*x(:,k) + a,2)];
    
end

figure(1);

semilogy([1:1:itr], grad_x)
title('Euclidean norm of gradient versus iteration index')
xlabel('Iteration index')
ylabel('Euclidean norm of gradient')
grid on

% x_err will save the distance from x_k to x^* at each iteration, in
% Euclidean norm

x_err = [];

for k=1:itr
   
    x_err = [x_err; norm(x(:,k)-x_opt,2)];
    
end

figure(2);

semilogy([1:1:itr], x_err)
title('Euclidean distance to the optimal solution versus iteration index')
xlabel('Iteration index')
ylabel('Distance to the optimal solution')
grid on
