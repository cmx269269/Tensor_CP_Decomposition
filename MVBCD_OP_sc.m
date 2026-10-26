function [F,count,fitchange,omega1,T,nonzerocols,normC] = MVBCD_OP_sc(XIJK,R,epsilon,alpha,itermax)
%MBCD矩阵型算法
% 构造缺失张量T
% load SD32_omega05.mat;
% load baizhu_omega05.mat;

omega1 = omega;
T = double(XIJK).*omega1;                 % 计算Hadamard积


m=size(T);
N=ndims(T);

G1=cell(N,1);
X=cell(N,1);
A=cell(N,1);
O=cell(N,1);
normX = norm(tensor(T));



A_in = cpd_gevd(double(T),R);
for i=1:N     % 随机生成初始组分以及将初始缺失张量T和缺失张量对应的张量omega进行模展开
    X{i}=tenmat(T,i).data;
    O{i}=tenmat(omega1,i).data;
    G1{i}=100*ones(m(i),R);
end
for i = 1:N-1
    A{i} = A_in{i};
end
A{N} = A_in{N};
M = m(1)*m(2)*m(3);
if M < 5e5
    n = 1;
else
    n = ceil(M/5e5)+11;
end
delta2 = norm((tensor(ktensor({A{1:N}}))-tensor(T)));
lambda =n*sqrt(2*delta2* log(200*M));
disp(lambda);
 
F=A;
F_bar = F;
m1 = 5; % 内迭代步数
nonzerocols = [];
F1 = 1e2;
normF = zeros(1,N);
normC = [];
for count = 1:itermax   % 最外层循环，用于计算最优的成分
    lambda = max(lambda * 0.99,1e-6); 
       
        A1 = A;
        A2 = A;
    for i = 1:N                % 内层循环，用于计算每一个组分以及更新
        C={A1{1:i-1},A1{i+1:end}}';
        D=khatrirao(C,'r');
        M = (O{i}.*X{i})*D;
        if i~=N   % 计算前N-1个成分
            F{i} = sub_BCD_VM(M,F{i},D,O{i},R,m1);
            A1{i} = F{i}+alpha(i)*(F{i}-A1{i});
            normF(i) = norm(F{i}-F_bar{i},'fro')^2;

        else      % 计算最后一个成分
         
            F{i} = sub_BCD_VMC1(M,F{i},D,O{i},R,m1,lambda);
            A1{i} = F{i}+alpha(i)*(F{i}-A1{i});
            normF(i) = norm(F{i}-F_bar{i},'fro')^2;
            nonzerocols(count) = sum(any(F{i},1));
            disp(nonzerocols(count));
        end
    
    end
    F2 = F_ABC(X,O,R,F,D,N) + 0.5 * normX^2;
    if  (F2 + 0.5 * sum(normF) *(epsilon/m1) > F1)  % 下降性检测
        F = F_bar;  % 重新采用原来的F进行下一次计算
        for i = 1:N
            C={F{1:i-1},F{i+1:end}}';
            D=khatrirao(C,'r');
            M = (O{i}.*X{i})*D;
            if i~=N   % 计算前N-1个成分
               F{i} = sub_BCD_VM(M,F{i},D,O{i},R,m1);
               A2{i} = F{i}+alpha(i)*(F{i}-A2{i});
            else
               F{i} = sub_BCD_VMC1(M,F{i},D,O{i},R,m1,lambda);
               A2{i} = F{i}+alpha(i)*(F{i}-A2{i});
            end

        end
        A = A2;
    else
        A = A1;
    end

    F1 = F2;
    F_bar = F;
    % 判断是否达到收敛条件
    fitchange = max(max(norm(G1{1}-F{1},'fro'),norm(G1{2}-F{2},'fro')),norm(G1{3}-F{3},'fro'));
    X_new = ktensor({F{1},F{2},F{3}});
    tensor_diff = norm(tensor(double(X_new) - double(XIJK)))/norm(tensor(double(XIJK)));
    disp(fitchange);
    crit = (fitchange < epsilon);
    if crit; nstall = nstall+1; else, nstall = 0; end
        if nstall >= 2
            break; 
        end


    G1 = F;
    disp(count);
    


end

end