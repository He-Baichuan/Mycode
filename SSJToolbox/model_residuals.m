function res = model_residuals(X, E, T, varargin)
% model_residuals  计算模型的残差
%
%   res = model_residuals(X, E, T, alpha, beta, gamma, delta, rho)
%   其中：
%     X     - (5T×1) 状态向量
%     E     - (T×1) 冲击向量
%     T     - 标本期数
%   可选参数（必须按顺序提供）：
%     alpha, beta, gamma, delta, rho

  % ----- 1. 拆包输入参数 -----
  if numel(varargin) < 5
    error('需要提供 alpha, beta, gamma, delta, rho 五个参数');
  end
  alpha = varargin{1};
  beta  = varargin{2};
  gamma = varargin{3};
  delta = varargin{4};
  rho   = varargin{5};

  % ----- 2. 预分配 -----
  res = zeros(5 * T, 1);

  % ----- 3. 逐期计算残差 -----
  for t = 1:T
    % 当前期
    Z = X(5*(t-1)+1);
    R = X(5*(t-1)+2);
    K = X(5*(t-1)+3);
    Y = X(5*(t-1)+4);
    C = X(5*(t-1)+5);

    % 滞后期
    if t == 1
      K_lag = K;
      Z_lag = Z;
    else
      K_lag = X(5*(t-2)+3);
      Z_lag = X(5*(t-2)+1);
    end

    % 前导期
    if t == T
      R_fwd = R;
      C_fwd = C;
    else
      R_fwd = X(5*t+2);
      C_fwd = X(5*t+5);
    end

    % 冲击
    epsilon = E(t);

    % 5 个方程的残差
    idx = 5*(t-1);
    res(idx+1) = -C^(-gamma) + beta * R_fwd * C_fwd^(-gamma);  % Euler
    res(idx+2) = -R + alpha * Z * K_lag^(alpha - 1) + 1 - delta;  % 利率
    res(idx+3) = -K + (1 - delta) * K_lag + Y - C;               % 资源约束
    res(idx+4) = -Y + Z * K_lag^alpha;                           % 生产函数
    res(idx+5) = -log(Z) + rho * log(Z_lag) + epsilon;          % TFP 冲击
  end
end
