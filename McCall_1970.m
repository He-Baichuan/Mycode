% 这是一个复现 McCall(1970) 的代码


n = 50;
alpha = 200;
beta = 100;
rng('default')
% 1. 支持集 (0 到 n)
k = 0:n;

% 2. 计算 Beta-Binomial 的 PMF
% P(X=k) = C(n,k) * B(k+α, n-k+β) / B(α, β)
% 为避免数值溢出，建议用对数计算后取指数
log_p = gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1) ...
      + betaln(k+alpha, n-k+beta) - betaln(alpha, beta);
p = exp(log_p)';

% 验证概率和为 1
%sum(p)  % 应接近 1.0

% 3. 线性间隔的工资向量 (10.0 到 60.0，共 n+1 = 51 个点)
w = linspace(10.0, 60.0, n+1)';

%%%%%%% 绘图 %%%%%%%
figure(1)
plot(w,p,'LineStyle','-','LineWidth',1.5)
xlabel('wage')
ylabel('Probability')
% 计算 样本均值和方差
E_w = dot(p,w)
Var_w = dot(p,w.^2)-E_w.^2
v_plot = zeros(n+1,7);
v_plot(:,1) = w/(1-beta);
for t = 2:7
    v_plot(:,t) = TV(v_plot(:,t-1),w,p);
end
figure(2)
plot(v_plot(:,2:7),'LineWidth',1.5)
legend('迭代次数1','迭代次数2','迭代次数3','迭代次数4','迭代次数5','迭代次数6')
beta = 0.99;
c = 25;
tol = 1e-8;
maxit = 500;
iter = 2;
crit = 1;
polcy_iter = [];
v_iter = [];
polcy_iter(:,1) = ones(n+1,1);
v_iter(:,1) = zeros(n+1,1);
coeff_matrix_accept = diag((1-beta)*ones(n+1,1));
coeff_matrix_search = -beta*repmat(p',[n+1,1])+eye(n+1);
tic
%%% 基于Howard 政策迭代 %%%
while iter<maxit && crit>tol
    %%% 政策评估 %%%
    index = (polcy_iter(:,iter-1)-1);
    A = index.*coeff_matrix_search+...
        (1-index).*coeff_matrix_accept;
    b = index.*(c*ones(n+1,1))+ (1-index).*w;
    v_iter(:,iter) = A\b;
    crit = max(abs(v_iter(:,iter)-v_iter(:,iter-1)),[],'all');
    % if ~any(polcy_iter(:,iter) ~= polcy_iter(:,iter-1))
    %     break;
    % end
    %%% 政策更新 %%%
    val_accept = w/(1-beta);
    val_search = c+beta*dot(p,v_iter(:,iter));
    polcy_iter(:,iter) = 1+(val_accept<val_search);
    iter = iter + 1;
end

w_R = (1-beta)*(c+beta*dot(p,v_iter(:,end)))
toc
figure(3)
plot(polcy_iter)
% 参数
% beta = 0.99; c = 25; n = 50; tol = 1e-6; maxit = 500;

% Beta-Binomial 概率（需先计算，或用你已有的 p）
% w = linspace(10, 60, n+1)';

% % 值迭代
% v = w / (1-beta);           % 初始值：全部接受
% error = inf; iter = 0;
% tic
% while error > tol && iter < maxit
%     % 贝尔曼算子 T(v)
%     accept_val = w / (1-beta);
%     reject_val = c + beta * dot(p, v);      % 标量
%     v_next = max(accept_val, reject_val);   % 逐点取max
% 
%     error = norm(v_next - v, inf);
%     v = v_next;
%     iter = iter + 1;
% end
% toc
% % 保留工资
% w_R = (1-beta) * (c + beta * dot(p, v))
% cdf = cumsum(p);
% H = 1 - interp1(w, cdf, w_R, 'spline');
% E_duration = (1/H)
dd = compute_mean_stopping_time(w_R, w, n, alpha, 100, 10000)



辅助函数

function [v_new] = TV(v_old,w,p)
beta = 0.99;
c = 25;
v_new = max(w/(1-beta),c+beta*dot(p,v_old));
end

function t = compute_stopping_time(w_bar, w, n, alpha, beta_param, seed)
    % 设置随机种子（Julia: Random.seed!）
    rng(seed);
    
    % 验证：保留工资不超过最大工资，且分布支撑匹配
    assert(w_bar <= w(end), '保留工资不能超过最大工资 w(end)');
    
    t = 1;
    while true
        % Beta-Binomial 两阶段抽样
        p_draw = betarnd(alpha, beta_param);   % 抽成功概率
        k = binornd(n, p_draw);                % 抽档位，取值 0,1,...,n
        
        % 关键：MATLAB 索引 1-based，k=0 对应 w(1)
        w_val = w(k + 1);
        
        if w_val >= w_bar
            return;   % MATLAB 中 return 直接退出并返回 t
        else
            t = t + 1;
        end
    end
end

function mean_time = compute_mean_stopping_time(w_bar, w, n, alpha, beta_param, num_reps)
    % 默认重复 10000 次（Julia 默认参数）
    if nargin < 6
        num_reps = 10000;
    end
    
    times = zeros(num_reps, 1);
    for i = 1:num_reps
        % 每次用不同的 seed（对应 Julia 的 seed = i）
        times(i) = compute_stopping_time(w_bar, w, n, alpha, beta_param, i);
    end
    
    mean_time = mean(times);
end
