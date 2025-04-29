clear; 
% clc;
tic;
%% 1. 参数与稳态计算（注意要先算出 K_ss、Y_ss、C_ss）
alpha = 1/3;
beta = 0.99;
theta = 1.7386;
delta = 0.02;
rho = 0.979;
T     = 400;

% 1.1 稳态值
nbar = 1/3;
kbar = nbar*(alpha/(1/beta-1+delta))^(1/(1-alpha));
Rkbar = 1/beta-1+delta;
wbar = (1-alpha)*(kbar/nbar)^(alpha);
nbar = (1+theta*((kbar/nbar)^(alpha)-delta*(kbar/nbar))/((1-alpha)*(kbar/nbar)^(alpha)))^(-1);
ybar = kbar^(alpha)*nbar^(1-alpha);
cbar = ybar-delta*kbar;
ibar = ybar-cbar;
rbar = 1/beta-1;
abar = 1;
par = struct('alpha',alpha,'beta',beta,'theta',theta,'delta',delta,'rho',rho);
%% 2. 构造 varFields（用 containers.Map 而非 struct）
vars.endo ={'a','r','k','y','c','Rk','w','i','n'};
vars.exog = {'epsilon'};
%% 3. 构造 steadystate 结构体
steadystate.initial  = [abar;rbar;kbar;ybar;cbar;Rkbar;wbar;ibar;nbar];
steadystate.terminal = steadystate.initial;
steadystate.exog     = 0;

%% 4. 实例化 ModelEnv
m = ModelUtils.ModelEnv(par,vars, steadystate, T);
[Xss, Ess] = ModelUtils.longsteadystate(m);
Ess(1) = 0.009; 
IRFmat = ModelUtils.linearIRFs(@RBC2_resid, m, Xss, Ess);
X_resp = Xss + IRFmat * Ess;
ModelUtils.plotIRFs(X_resp, m, T, 1, 0);
toc;










