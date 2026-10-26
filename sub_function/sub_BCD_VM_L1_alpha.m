function [X] = sub_BCD_VM_L1_alpha(T,X,Q,O,s,lambda) 
% X表示计算的α
% O表示缺失张量的投影
% T 表示缺失张量
% s 表示大于1的常数
% eta = max(norm(Q,'fro')^2,1);
eta = max(norm(Q*Q','fro'),1);
Ovec = vec(tensor(O));
Tvec = vec(T);
gradX = ((Ovec').*(X*Q)- Tvec')*Q';
X = X - gradX/(s*eta);
%% 软阈值算子
X = sign(X).*max(abs(X)-lambda/(s*eta),0);
    



end