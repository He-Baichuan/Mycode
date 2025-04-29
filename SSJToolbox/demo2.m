clear; 
% clc;
tic;
%% 1. 参数与稳态计算（注意要先算出 K_ss、Y_ss、C_ss）
alpha = 0.36; 
beta  = 0.96; 
gamma = 2.0; 
delta = 0.08; 
rho   = 0.9;
T     = 180;
pars = struct('alpha',alpha,'beta',beta,'gamma',gamma,'delta',delta,'rho',rho);
% 1.1 稳态值
Z_ss = 1.0;
R_ss = 1 / beta;
K_ss = ((R_ss - 1 + delta)/alpha)^(1/(alpha - 1));
Y_ss = K_ss^alpha;
C_ss = Y_ss - delta * K_ss;
  
%% 2. 构造 varFields（用 containers.Map 而非 struct）

vars.endo = {'Z','R','K','Y','C'};
vars.exog = {'epsilon'};

%% 3. 构造 steadystate 结构体
steadystate.initial  = [Z_ss; R_ss; K_ss; Y_ss; C_ss];
steadystate.terminal = steadystate.initial;
steadystate.exog     = 0;

%% 4. 实例化 ModelEnv
m = ModelUtils.ModelEnv( pars, vars, steadystate, T);

%% 5. 测试 steady-state 取值
% 比如取 C 的 steady state
% disp(ModelUtils.steadystatevalue('C', m, 'endog', 'initial'));

ModelUtils.checksteadystate(m, @rbc_resid);
% % 3.1 构造稳态路径
[Xss, Ess] = ModelUtils.longsteadystate(m);
% 
% % 3.2 用类自带的 finite‐difference Jacobian
% %     linearIRFs 内部会调用 get_fX 和 get_fE
Ess(1) = 0.01; 
IRFmat = ModelUtils.linearIRFs(@rbc_resid, m, Xss, Ess);
% FX = ModelUtils.get_fX(@rbc_resid, m, Xss, Ess);

% 
% % 3.4 得到响应路径：X_resp = Xss + IRF * shock
X_resp = Xss + IRFmat * Ess;
% % 用 wide 拆分

ModelUtils.plotIRFs(X_resp, m, T, 1, 0);

toc;









