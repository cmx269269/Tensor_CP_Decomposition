clc;clear;close all;
% load SD32.mat;
% load baizhu.mat
xx11 = X;
C_real=C;
itermax = 5000;
xx11=tensor(xx11);
X1=tensor(xx11);
epsilon=1e-6;   
R=30;
alpha=[0.3,0.5,0.8];
t = [];
err = [];
i1=0;
j1=0;
k1=0;
            for kk = 1:1
                tic
                [D1,count,fitchange,omega,T,nonzerocols,normC]=MVBCD_OP_sc(xx11,R,epsilon,alpha,itermax);  % 计算前N个因子矩阵
                t(kk)=toc;
                sol1=double(ktensor(D1)); 
                X_new = double(ktensor(D1));
                tensor_diff = norm(tensor(double(X_new) - double(xx11)))/norm(tensor(double(xx11)));
                C1=D1{3}; % 计算浓度矩阵
                A1 = D1{1};
                B1 = D1{2};
                save('mydata','A1','B1','C1','tensor_diff');
                [RE_ER,eerr,C_last] = huishoulv2(C_real,abs(C1));
            end
              %   % RMSEP = sum(eer)/10; 
              % t1=sum(t)/10;    % 计算时间平均值
              % ER=sum(eer)/10;   % 计算误差平均值

% t2 = tenmat(t1,1);
% [m,n] = size(double(t2));
% x = 1:m;
% plot(x,t3,'-*y',x,d2(:,1),'-*r',x,d2(:,2),'-*g',x,d2(:,3),'-*c',x,d2(:,4),'-*b',x,d2(:,5),'-*m',x,d2(:,6),'-*k')
% ylabel('Time(s)');
% xlabel('\beta_1');
% legend('\beta_2 = 0,\beta_3 = 0','\beta_2 = 0.4,\beta_3 = 0.7','\beta_2 = 0.4,\beta_3 = 0.8','\beta_2 = 0.5,\beta_3 = 0.7','\beta_2 = 0.5,\beta_3 = 0.8','\beta_2 = 0.6,\beta_3 = 0.7' ...
%     ,'\beta_2 = 0.6,\beta_3 = 0.8');
[m,~] = size(C1);
x = 1:m;
% PNSR = 10 * log10(a*b*c*norm(X1,'inf')/norm(tensor(D1)-X1));
if R == 3
% figure(1);
% plot(x,C_last(:,1),'-or',x,C_last(:,2),'-hg');
% legend('Clast-1','Clast-2');
% figure(2);
% plot(x,C_real(:,1),'-or',x,C_real(:,2),'-hg');
% legend('C1-1','C1-2');
figure(3);
C2 = abs(C1);
plot(x,C2(:,1),'-or',x,C2(:,2),'-hg',x,C2(:,3),'-*c');
xlabel('Sample number');
ylabel('Relative concentration');
legend('com 1','com 2','com 3');
elseif R == 5
% figure(1);
% plot(x,C_last(:,1),'-or',x,C_last(:,2),'-hg');
% legend('Clast-1','Clast-2');
% 
% figure(2);
% plot(x,C_real(:,1),'-or',x,C_real(:,2),'-hg');
% legend('C1-1','C1-2');

figure(3);
C2 = abs(C1);
plot(x,C2(:,1),'-or',x,C2(:,2),'-hg',x,C2(:,3),'-*c',x,C2(:,4),'-db',x,C2(:,5),'-c');
xlabel('Sample number','FontSize',16);
ylabel('Relative concentration','FontSize',16);
legend('com 1','com 2','com 3','com 4','com 5');
elseif R == 7
% figure(1);
% plot(x,C_last(:,1),'-or',x,C_last(:,2),'-hg');
% legend('Clast-1','Clast-2');
% figure(2);
% plot(x,C_real(:,1),'-or',x,C_real(:,2),'-hg');
% legend('C1-1','C1-2');
figure(3);
C2 = abs(C1);
plot(x,C2(:,1),'-or',x,C2(:,2),'-hg',x,C2(:,3),'-*c',x,C2(:,4),'-db',x,C2(:,5),'-c',x,C2(:,6),'-ok',x,C2(:,7),'-ob');
xlabel('Sample number');
ylabel('Relative concentration');
legend('com 1','com 2','com 3','com 4','com 5','com 6','com 7');
elseif R==9
figure(3);
C2 = abs(C1);
plot(x,C2(:,1),'-or',x,C2(:,2),'-hg',x,C2(:,3),'-*c',x,C2(:,4),'-db',x,C2(:,5),'-c',x,C2(:,6),'-ok',x,C2(:,7),'-ob',x,C2(:,8),'-b',x,C2(:,9),'-og');
xlabel('Sample number');
ylabel('Relative concentration');
legend('com 1','com 2','com 3','com 4','com 5','com 6','com 7','com 8','com 9');


else
    disp("无");

end


% ,x,C2(:,4),'-db',x,C2(:,5),'-c',x,C2(:,6),'-ok',x,C2(:,7),'-ob'
% figure(1);
% plot(x,C_last(:,1),'-or',x,C_last(:,2),'-hg',x,C_last(:,3),'-*c');
% legend('Clast-1','Clast-2','Clast-3');
% % figure(1);
% % plot(C1,'-o');
% figure(2);
% plot(x,C_real(:,1),'-or',x,C_real(:,2),'-hg',x,C_real(:,3),'-*c');
% legend('C1-1','C1-2','C1-3');
% % figure(2);
% % plot(x,C1(:,1),'-or',x,C1(:,2),'-hg',x,C1(:,3),'-*c',x,C1(:,4),'-db',x,C1(:,5),'-c',x,C1(:,6),'-ok');
% % legend('C1-1','C1-2','C1-3','C1-4','C1-5','C1-6');

% plot(x,C1(:,1),'-or',x,C1(:,2),'-hg',x,C1(:,3),'-*c',x,C1(:,4),'-db');
% legend('C1-1','C1-2','C1-3','C1-4');