clc;
clear;
%% 参数设定
r = 0.01;
y1 = 1;
y2 = 2;
beta = 0.98;
a_lower = 0;%borrowing limit
P = [0.3,0.7;0.3,0.7];%s的概率转移矩阵
N = 200;%gridpoint数目
epsilon = 0.0001;%tolerable error

%% 迭代过程
u = linspace(0,1.5,N);
a = exp(exp(u)-1)-1+a_lower;%格点
A = [a;a];
Y = [ones(1,N)*y1;ones(1,N)*y2];
C = ((1+r)*A+Y)/2;%猜测的policy funtion，作为初始值。注意不要破坏budget constraint

delta = 1;
while delta >= epsilon
    Vp = (1+r)./C;
    C_endo = 1./(beta*P*Vp);
    F_ap = (A+C_endo-Y)/(1+r);
    Ap = [interp1(F_ap(1,:),a,a,'linear','extrap');interp1(F_ap(2,:),a,a,'linear','extrap')];
    Ap = max(Ap,a_lower);
    C_new = (1+r)*A+Y-Ap;
    delta = max(abs(C_new-C),[],'all');
    C = C_new;
end
figure ;
plot(a,Ap(1,:),a,Ap(2,:),a,a);
figure;
plot(a,C(1,:),a,C(2,:));