clear;
clc;

B = [12.6695,0,-1.2353,0,0;0,1,0,0,0;0,-1,0.36,0,0;0,0,1,0,0;0,0,0,1,-0.03475];
A = [12.353,0,0,-0.9186,0;0,0.95,0,0,0;0.36,0,0,-0.64,0;1,0,0,0,1;0,0,0,1,0];
G = [0;1;0;0;0];
m = 2;%非前定变量的个数是m
[C, D, N, L, eigval] = BKsolve(A, B, G, m);
% [ZZ,TT,SS,alpha,beta,QQ] =schurg(A,B);
% eigvalue = alpha./beta;
% Zp = ZZ';
% [a,b] = size(ZZ);
% 

% N = inv(Zp(a-m+1:a,a-m+1:a))*Zp(a-m+1:a,1:a-m);
% BBN = inv(B(1:a-m,1:a-m)-B(1:a-m,a-m+1:a)*N);
% C = BBN*(A(1:a-m,1:a-m)-A(1:a-m,a-m+1:a)*N);
% L1 = inv(Zp(a-m+1:a,a-m+1:a))*inv(SS(a-m+1:a,a-m+1:a));
% L2 = QQ(a-m+1:a,1:a-m)*G(1:a-m,1)+QQ(a-m+1:a,a-m+1:a)*G(a-m+1:a,1);
% L = L1*L2;
% D = BBN*(G(1:a-m,1)-A(1:a-m,a-m+1:a)*L);


T = 120;
tt = (1:T);
shock = zeros(1,T);
shock(1) = 1;
Cs = size(C,1);
Ns = size(N,1);
x_irf = zeros(Cs,T);
y_irf = zeros(Ns,T);

for i=1:T
    if i== 1
        x_irf(:,i) = D*shock(i);
        y_irf(:,i) = -L*shock(i);
    else
        x_irf(:,i) = C*x_irf(:,i-1) + D*shock(:,i);
        y_irf(:,i) = -N*x_irf(:,i-1)-L*shock(:,i);
    end
end
xy_irf = [x_irf;y_irf];
figure
plot(tt,xy_irf,tt,zeros(1,T));



