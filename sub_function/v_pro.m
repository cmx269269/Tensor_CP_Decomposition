function [x_hat] = v_pro(x)     % 单位约束的邻近算子
if norm(x) ~= 0
    x_hat = x/norm(x);
else
    x(1) = 1;
    x_hat = x;
end

end