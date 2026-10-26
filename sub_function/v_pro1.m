function [x_hat] = v_pro1(x,lambda)     % 20范数的邻近算子
epsilon1 = 1e-4;
L = length(x);
if norm(x) > sqrt(2*lambda)
    x_hat = x;
elseif norm(x) == sqrt(2*lambda)
    x_hat = zeros(L,1);
else
    x_hat = zeros(L,1);
end

end