clear;
clc;
data = readmatrix('/Users/mac/Desktop/dataset/中国季度数据.xlsx','Sheet', ...
    'Data','Range','D2:D97');%读取数据
Y = data(:,1);
Y = log(Y);
T = size(Y,1);
T = T(1);
t = (1:T)';
tt = (2000:.25:2023.75)';
x = [ones(T,1),t,t.^2];
coe = inv(x'*x)*(x'*Y);
a = coe(1);
b = coe(2);
c = coe(3);
y_hat = a*x(:,1)+b*x(:,2)+c*x(:,3);
y_c = Y-y_hat;
yp = 100 * (exp(y_c) - 1);
% yp = 100 * y_c;
ySigma = std(y_c);
% t = 1999.75+t;

figure('Position', [100, 100, 1200, 400]);
subplot(1,2,1);
plot(tt,y_hat,'-k',tt,Y,'-.k','LineWidth',1);
xlabel('年份');
% ylabel('百分比');
title(' 中国真实人均GDP的趋势与周期');
xlim([2000, tt(end)]);

subplot(1,2,2)
plot(tt,yp,'-k','LineWidth',1);
xlabel('年份');
ylabel('趋势偏离的百分比');
title(' 中国真实人均GDP的趋势与周期');
xlim([2000, tt(end)]);


