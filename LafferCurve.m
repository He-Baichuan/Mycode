clear;
clc;

%参数校准
par.alpha = 0.38;
par.delta = 0.07;
par.nu0 = 3.732;
par.nu1 = 1;
par.sigma = 2;
par.Lbar = 0.25;
par.tau_l = 0.28;
par.tau_k = 0.36;
par.tau_c = 0.05;
par.r_B = 0.04;
par.psi = 1.02;
par.BtoY = 0.63;
par.GtoY = 0.18;

YK = par.r_B/(par.alpha * (1 - par.tau_k))+ par.delta/par.alpha;
YL = YK^(-par.alpha/(1-par.alpha));
Y = YL * par.Lbar;
K = Y/YK;
w = (1-par.alpha)*YL;
RHS = 1/(par.sigma * par.nu0 * par.Lbar^(1+1/par.nu1))+1 -1/par.sigma;
CY = (1-par.tau_l)/(1+par.tau_k)*(1-par.alpha)/(1+1/par.nu1)*RHS;
C = CY*Y;
I = (par.psi-1+par.delta)*K;
G = Y- I - C;
x0 = [K;C;Y;par.Lbar;I];
epsilon = 1e-1;
par.tau_l = par.tau_l-epsilon;
obj = @(x) findsteadystate(x,par);
x1 = fsolve(obj,x0);
w = (1-par.alpha)*x1(3)/x1(4);
taul = par.tau_l;
tt1 = w*taul*x1(4) + par.tau_k*(par.r_B/(1-par.tau_k)+par.delta )*x1(1) + par.tau_c*x1(2);

par.tau_l = par.tau_l+2*epsilon;
x1 = fsolve(obj,x0);
w = (1-par.alpha)*x1(3)/x1(4);
taul = par.tau_l;
tt2 = w*taul*x1(4) + par.tau_k*(par.r_B/(1-par.tau_k)+par.delta )*x1(1) + par.tau_c*x1(2);
T = (tt2-tt1)/(2*epsilon);%自我融资率
n_tau = 99;
seq_tau_l = linspace(0.01,0.99,n_tau);
result = zeros(n_tau,3);
for i = 1:n_tau
    taul = seq_tau_l(i);
    result(i,1) = 100*taul;
    par.tau_l = taul;
    obj = @(x) findsteadystate(x,par);
    x1 = fsolve(obj,x0);
    w = (1-par.alpha)*x1(3)/x1(4);
    result(i,2) = w*taul*x1(4);
    result(i,3) = w*taul*x1(4) + par.tau_k*(par.r_B/(1-par.tau_k)+par.delta )*x1(1) + par.tau_c*x1(2);
    x0 = x1;
end
[val1, id1] = max(result(:,2));
[val2, id2] = max(result(:,3));
figure
plot(result(:,1),result(:,2),'linewidth',1.5,'color',[178/255,34/255,34/255] )
hold on
plot(result(:,1),result(:,3),'linewidth',1.5,'color',[8/255,62/255,118/255] )
xline(28,'LineWidth',1,'LineStyle','--')
plot(result(id1),val1,'kx', 'MarkerSize', 12,'linewidth',1.5)
plot(result(id2),val2,'kx', 'MarkerSize', 12,'linewidth',1.5)
legend('Labor Income Tax','All tax')
ylim([0,0.18])
xlim([0,100])

seq_tau_k = linspace(0.01,0.99,n_tau);
result2 = zeros(n_tau,3);
par.tau_l = 0.28;
for i = 1:n_tau
    tauk = seq_tau_k(i);
    result2(i,1) = 100*tauk;
    par.tau_k = tauk;
    obj = @(x) findsteadystate(x,par);
    x1 = fsolve(obj,x0);
    w = (1-par.alpha)*x1(3)/x1(4);
    r = (par.r_B/(1-par.tau_k)+par.delta );
    result2(i,2) = r*tauk*x1(1);
    result2(i,3) = w*par.tau_l*x1(4) +r*tauk*x1(1) + par.tau_c*x1(2);
    x0 = x1;
end
[val1, id1] = max(result2(:,2));
[val2, id2] = max(result2(:,3));
figure
plot(result2(:,1),result2(:,2),'linewidth',1.5,'color',[178/255,34/255,34/255] )
hold on
plot(result2(:,1),result2(:,3),'linewidth',1.5,'color',[8/255,62/255,118/255] )
xline(36,'LineWidth',1,'LineStyle','--')
plot(result(id1),val1,'kx', 'MarkerSize', 12,'linewidth',1.5)
plot(result(id2),val2,'kx', 'MarkerSize', 12,'linewidth',1.5)
legend('Capital Income Tax','All tax')
ylim([0,0.18])
xlim([0,100])


function [y] = findsteadystate(x,par)
Kss = x(1);
Css = x(2);
Yss = x(3);
Lss = x(4);
Iss = x(5);
Gss = par.GtoY*Yss;
y(1) = Yss - Css - Gss - Iss;
YK = par.r_B/(par.alpha * (1 - par.tau_k))+ par.delta/par.alpha;
YL = YK^(-par.alpha/(1-par.alpha));
Y = YL * Lss;

K = Y/YK;
% w = (1-par.alpha)*YL;
RHS = 1/(par.sigma * par.nu0 * Lss^(1+1/par.nu1))+1 -1/par.sigma;
CY = (1-par.tau_l)/(1+par.tau_k)*(1-par.alpha)/(1+1/par.nu1)*RHS;
y(2) = CY - Css/Yss;
I = (par.psi-1+par.delta)*K;
% G = Y- I - C;
y(3) = Yss - Y;
y(4) = K-Kss;
y(5) = I - Iss;
% y = [K;C;Y;L;I];
end
