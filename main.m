close all
clear
clc

% Load the data
xtr2001 = csvread('data_matlab/xtr2001.csv', 1,1); % dim 1: iso_o; dim 2: iso_d
tar2001 = csvread('data_matlab/tariff2001.csv', 1,1); % dim 1: iso_o; dim 2: iso_d

xtr = xtr2001;
tar = tar2001;

m.N = 36;
m.sigma = 3.059;
m.theta = m.sigma - 1;
%%%% 利用data计算的值
X = sum(xtr,1)';% X_N
lambda = xtr./repmat(X',[m.N,1]);
wL = sum(1./(1 + tar).*xtr,2);

m.X = X;
m.lambda = lambda;
m.wL = wL;
m.tar = tar;

ttar = m.tar;
htau = ones(m.N, m.N);
%%%% 计算exact hat algebra，以便校准数据（遭受冲击后），因为数据不一定与Z'相一致
tic
[hw, hP, hX, hlambda, hwel] = eqm_iter(ttar, htau, m);
toc

m.lambda = m.lambda.*hlambda; % Balanced trade share
m.X = m.X.*hX; % Balanced trade expenditure
m.wL = m.wL.*hw; % Balanced trade wage

% Test the trade balance
[hw, hP, hX, hlambda, hwel] = eqm_iter(ttar, htau, m);
%% Counterfactual: China eliminates all of its import tariffs

ttar = m.tar;
i = 9; % China
ttar(:,i) = zeros(m.N, 1);
htau = ones(m.N, m.N);

[hw, hP, hX, hlambda, hwel] = eqm_iter(ttar, htau, m);

dwel = (hwel-1)*100;
drw = (hw./hP - 1)*100;
%% Counterfactual: China's unilateral optimal import tariffs

i = 9; % China
options = optimset('Display', 'iter','MaxIter', 10000, 'MaxFunEval', 100000000000 ,...
        'TolX', 1e-16, 'TolFun', 1e-16, 'TolCon', 1e-10, 'Algorithm', 'sqp');
% Lower and upper bound
hX_lb = 0*ones(m.N, 1);
hP_lb = 0*ones(m.N, 1);
hw_lb = 0*ones(m.N, 1);
ttar_lb = -1*ones(m.N, 1);
ttar_lb(i) = 0;
lb = [hX_lb; hP_lb; hw_lb; ttar_lb];

hX_ub = inf*ones(m.N, 1);
hP_ub = inf*ones(m.N, 1);
hw_ub = inf*ones(m.N, 1);
ttar_ub = inf*ones(m.N, 1);
ttar_ub(i) = 0;
ub = [hX_ub; hP_ub; hw_ub; ttar_ub];

x0 = [ones(3*m.N, 1); m.tar(:,i)]; % Initial guesses
tic
[x,fval,exitflag] = fmincon(@(x)obj_optimalTariff(x, i, m),x0,[],[],[],[],lb,ub,@(x)constraint_optiTariff(x, i ,m), options);
toc

% hX = x(1: m.N);
% hP = x(m.N + 1: 2*m.N);
% hw = x(2*m.N + 1: 3*m.N);
ttar_cn = x(3*m.N + 1: 4*m.N);

ttar = m.tar;
ttar(:,i) = ttar_cn;
htau = ones(m.N, m.N);
[hw, hP, hX, hlambda, hwel] = eqm_iter(ttar, htau, m);
dwel_opt = (hwel-1)*100;

disp('ttar_cn=')
disp(ttar_cn')
disp('dwel_cn = ')
disp(dwel_opt(i))