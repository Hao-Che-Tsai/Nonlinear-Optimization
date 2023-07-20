load('HW5_data.mat');
mu = 1e0;

n = length(c);
x_0 = zeros(n,1);
x_star =  gaussseidel(Q+mu*c*c',b,x_0,1e-9,3e5);

error = norm((Q+mu*c*c')*x_star-b)
