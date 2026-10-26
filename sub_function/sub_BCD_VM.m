function [X] = sub_BCD_VM(M,X,D,O,n,m)       
% n表示需要计算的变量个数
% m表示内迭代步数
% M表示对应组分矩阵与其余矩阵的hadmard积相乘
% X表示需要计算的组分
% O表示投影的模展开
epsilon1 = 1e-6;
Y=D'*D;  % Y表示对应的L构成的系数矩阵
for k = 1:m    % 计算前N-1个成分m次
    for j = 1:n  %  用于计算前N-1个成分的R个列向量
        gradX = (O.*(X*D'))*D(:,j) - M(:,j);
        X(:,j) = X(:,j) - gradX/(Y(j,j)+epsilon1);
        X(:,j) = v_pro(X(:,j));
    end
end

end