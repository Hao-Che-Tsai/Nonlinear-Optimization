% Implementation of Gradient Descent with Goldstein rule for inexact line
% search to solve
% min_x, y f(x,y) = x^2 - 5xy + y^4 - 25x - 8y


%% Gradient Descent

itr = 1;
x = [0; 0;]; % initial point
grad_x = norm(obj_grad(x),2); % grad_x will save Euclidean norm of gradients
f_x = obj(x); % f_x will save objective function values
tau = 10^(-6);

while (norm(obj_grad(x(:,itr)),2) > tau) % Stopping criterion: Euclidean norm of gradient <= epsilon
    
    x_itr = x(:,itr);
    g = obj_grad(x_itr)';
    d = -g; % Gradient Descent
    
    
    % Goldstein rule below
    % The code below will search for the step size that satisfies the
    % Goldstein criterion.
    
    % Initial interval of the step size: we will keep reducing this
    % interval until we find a step size satisfying Goldstein rule.
    alpha_u = 1;
    alpha_l = 0;
    %    eta = 1.2;
    epsilon = 0.4;
    
    
    % Below, increase alpha_u, until the lower bound of Goldstein criterion
    % is satisfied.
    while (obj(x_itr + alpha_u*d) < (obj(x_itr) + (1-epsilon)*alpha_u*g'*d))
        alpha_u = alpha_u*2;
    end
    
    % Note that we will keep alpha_u such that it will always satisfy the
    % lower bound of Goldstein rule in every update
    
    pass = 0;
    
    while (pass == 0)
        if (obj(x_itr + alpha_u*d) > (obj(x_itr) + (epsilon)*alpha_u*g'*d))
            
            alpha_tp = (alpha_u+alpha_l)/2;
            
            if (obj(x_itr + alpha_tp*d) <= (obj(x_itr) + (epsilon)*alpha_tp*g'*d)) && (obj(x_itr + alpha_tp*d) >= (obj(x_itr) + (1-epsilon)*alpha_tp*g'*d))
                alpha = alpha_tp;
                pass = 1;
            elseif (obj(x_itr + alpha_tp*d) > (obj(x_itr) + (epsilon)*alpha_tp*g'*d)) && (obj(x_itr + alpha_tp*d) >= (obj(x_itr) + (1-epsilon)*alpha_tp*g'*d))
                alpha_u = alpha_tp;
            else
                alpha_l = alpha_tp;
            end
        else
            alpha = alpha_u;
            pass = 1;
        end
    end
    
    alpha
    
    x_new = x_itr + alpha*d; % Gradient Descent step
    x = [x x_new];
    itr = itr + 1;
    
    grad_x = [grad_x, norm(obj_grad(x_new),2)];
    f_x = [f_x, obj(x_new)];
    
    
    itr
    
end

figure();

semilogy([1:1:itr], grad_x)
title('Euclidean norm of gradient versus iteration index')
xlabel('Iteration index')
ylabel('Euclidean norm of gradient')
grid on


figure();

plot([1:1:itr], f_x)
title('Objective function value versus iteration index')
xlabel('Iteration index')
ylabel('Objective function value')
grid on


x_opt = [20; 3] % Numerically, you can find that this would be the local minimum point to which the algorithm converges to.

% x_err will save the distance from x_k to x^* at each iteration, in
% Euclidean norm

x_err = [];

for k=1:itr
    
    x_err = [x_err; norm(x(:,k)-x_opt,2)];
    
end

figure();

semilogy([1:1:itr], x_err)
title('Euclidean distance to the relative minimum point versus iteration index')
xlabel('Iteration index')
ylabel('Distance to the optimal solution')
grid on

% Output the Hessian at the final solution point
obj_Hessian(x(:,itr))