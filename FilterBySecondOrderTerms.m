clear;
clc;
% 读取Excel文件（数值部分）
data = readmatrix('/Users/mac/Desktop/dataset/中国宏观经济季度数据.xlsx','Range','C2:K129');  % 忽略文本和原始数据
Y = data(:,1);
C = data(:,2);
% Y = log(Y);
T = size(Y);
T2 = size(C,1);
T = T(1);
t = (1:T)';
t2 = (1:T2)';
tt = (1992:0.25:2023.75);
tt2 = (1995:0.25:2017.75);

[y2_hat,y2_c] = hpfilter(Y,1600);
y2_p = 100*(exp(y2_c)-1);
[C2_hat,C2_c] = hpfilter(C,1600);

x = [ones(T,1),t,t.^2,t.^3];
x2 =[ones(T2,1),t2,t2.^2,t2.^3];
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
ySigma2 = std(y2_c);
coe2 = inv(x2'*x2)*(x2'*C);
a2 = coe2(1);
b2 = coe2(2);
c2 = coe2(3);
d2 = coe2(4);
C_hat = a2*x2(:,1)+b2*x2(:,2)+c2*x2(:,3)+d2*x2(:,4);
C_c = C-C_hat;
C_p = 100*(exp(C_c)-1);
C2_p = 100*(exp(C2_c)-1);
CSigma = std(C_c);
CSigma2 = std(C2_c);
Cacf = autocorr(C_c);
ration = CSigma/ySigma;


figure('Position', [100, 100, 1200, 400]);
subplot(1,2,1);
plot(tt,C_hat,'r',tt,y_hat,'-k','LineWidth',1);
xlabel('年份');
ylabel('趋势的对数值');
title(' 中国真实人均GDP的趋势与周期');
xlim([1992, tt(end)]);

subplot(1,2,2)
plot(tt,zeros(1,T),tt,C_p,'-r',tt,y_p,'-k','LineWidth',1);
xlabel('年份');
ylabel('趋势偏离的百分比');
title(' 中国真实人均GDP的趋势与周期');
xlim([1992, tt(end)]);
