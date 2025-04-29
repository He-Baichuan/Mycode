clear;
clc;
data1 = readmatrix('D:\data\债务数据.xlsx','Sheet','Sheet2','Range','B28:B45');
Y = data1;
Y = log(Y);
T = size(Y,1);
t = (1:T)';
tt = (2006:1:2023);
x = [ones(T,1),t,t.^2,t.^3];
coe = inv(x'*x)*(x'*Y);
a = coe(1);
b = coe(2);
c = coe(3);
d = coe(4);
y_hat = a*x(:,1)+b*x(:,2)+c*x(:,3)+d*x(:,4);
y_c = Y-y_hat;
y_p = 100*(exp(y_c)-1);
ySigma = std(y_c);
yacf = autocorr(y_c);

figure('Position', [100, 100, 1200, 400]);
subplot(1,2,1);
plot(tt,y_hat,'-k','LineWidth',1);
xlabel('年份');
ylabel('趋势的对数值');
title(' 中央财政债务余额的趋势与周期');
xlim([2006, tt(end)]);

subplot(1,2,2)

plot(tt,zeros(1,T),tt,y_p,'-k','LineWidth',1);
xlabel('年份');
ylabel('趋势偏离的百分比');
title(' 中央财政债务余额的趋势与周期');
xlim([2006, tt(end)]);

