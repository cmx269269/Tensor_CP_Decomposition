function [X] = sub_BCD_VMC1(M,X,D,O,n,m,lambda)
% n表示需要计算的变量个数
% m表示内迭代步数
% M表示对应组分矩阵与其余矩阵的hadmard积相乘
% X表示需要计算的组分
% Y表示对应的L构成的系数矩阵
% O表示投影的模展开
Y=D'*D;
for k = 1:m       % 计算最后一个成分m次
    for j = 1:n   % 用于计算最后一个成分的R个列向量
        gradX = (O.*(X*D'))*D(:,j) - M(:,j);
        X(:,j) = X(:,j) - gradX/Y(j,j);
        X(:,j) = v_pro1(X(:,j),lambda);
    end
end

end
