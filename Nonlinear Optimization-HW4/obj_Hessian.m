function F = obj_Hessian(x)
    
    F = zeros(2,2);
    F(1,1) = 2;
    F(1,2) = -5;
    F(2,1) = -5;
    F(2,2) = 12*x(2)^2;
end