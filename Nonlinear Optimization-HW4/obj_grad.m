function f_grad = obj_grad(x)
     
    f_grad(1) = 2*x(1) - 5*x(2) - 25;
    f_grad(2) = -5*x(1) + 4*x(2)^3 -8;
end