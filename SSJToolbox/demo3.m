clear;
clc;
%% 参数设定
alpha = 0.4;
beta = 0.98;
gamma = 2;
delta = 0.02;
rho = 0.95;
par =struct('alpha',alpha,'beta',beta,'gamma',gamma,'delta',delta,'rho',rho);



grid = HAGrids();
% disp(['资产网格点数：', num2str(grid.na)]);
% disp(['冲击状态数：',   num2str(grid.ne)]);
% disp('平稳分布 e:');    disp(grid.e');
% disp('转移矩阵 Pi_e:'); disp(grid.Pi_e);


obj2 = @(K) beta * getRonly(par,K) - 1; 
[Kss2, fval, exitflag] = fsolve(obj2, 48);
[W,R] = getWR(par,Kss2);
c0 = (R-1) * grid.a + W*grid.e';
ss = containers.Map();   % 创建一个空 Map
ss('c') = c0;            % 用字符串 'c' 作为 key，把 c0 赋给它

[Va, g, c] = SolveEGM(ss('c'),[W,R],par,grid);
% 假设 grid.a 是 na×1 向量，c 是 na×2 矩阵

figure;
plot(grid.a, c(:,1), 'LineWidth', 1.5);   % 第一列，对应低 endowment
hold on;
plot(grid.a, c(:,2), 'LineWidth', 1.5);   % 第二列，对应高 endowment
hold off;

title('Consumption Policy Rules');
legend('Low endowment', 'High endowment', 'Location', 'best');
xlabel('Assets');
ylabel('Consumption');

aBar = 100;  % 只画前 100 点
figure;
% 低 endowment 的储蓄决策（红实线）
plot(grid.a(1:aBar), g(1:aBar,1), 'r-', 'LineWidth', 1.5);
hold on;
% 高 endowment 的储蓄决策（蓝实线）
plot(grid.a(1:aBar), g(1:aBar,2), 'b-', 'LineWidth', 1.5);
% 45° 参考线（黑虚线）
plot(grid.a(1:aBar), grid.a(1:aBar), 'k--', 'LineWidth', 1.5);
hold off;

title('Savings Policy Rules');
legend('Low endowment', 'High endowment', '45-degree', 'Location', 'best');
xlabel('Assets');
ylabel('Savings');






%% 函数块
function val = uPrime(par,c)
    val = c.^(-par.gamma);
end

function val = uPrimeInv(par,up)
    val = up.^(-1.0/par.gamma);
end

function [W, R] = getWR(par,K)
    alpha = par.alpha;
    delta = par.delta;
    Z = 1;
    Lbar = 1.0;
    W = (1-alpha)*Z.*(K /Lbar).^alpha;
    R = alpha*Z.*(K/Lbar).^(alpha-1) + 1 - delta;   
end

function r = getRonly(par, K)
    [~, r] = getWR(par, K);
end

function val = aggregator(g,c,D)
    c_vec = c(:);
    val =  kron(D',g)*c_vec;
end

function idx = lookup(x, x1)
% LOOKUP  找出区间索引：对每个 x1 值，返回 i 使得 x(i-1) < x1 <= x(i)
%   输入：
%     x  - 已排序的列向量，长度 N
%     x1 - 待查找的值向量
%   输出：
%     idx - 与 x1 同大小的索引向量，取值范围 [2, N]
%
% 算法：
%   对每个 x1(j)：
%     k = 最后一个满足 x(k) <= x1(j) 的位置
%     将 k 限制到 [1, N-1]
%     idx(j) = k + 1

    N   = numel(x);
    idx = zeros(size(x1));
    for j = 1:numel(x1)
        % 找到最后一个 x(k) <= x1(j)
        k = find(x <= x1(j), 1, 'last');
        if isempty(k)
            k = 1;
        end
        % 限制 k 在 [1, N-1]
        k = min(max(k, 1), N-1);
        idx(j) = k + 1;
    end
end

function [y1, idx] = interp_custom(x, y, x1, left, right)
% INTERP_CUSTOM  一维线性插值/外推 (省略字符版)
%
% 用法：
%   y1 = interp_custom(x, y, x1)            % 默认 left=0, right=0
%   y1 = interp_custom(x, y, x1, 1)          % left=1, right=0
%   y1 = interp_custom(x, y, x1, 1, 1)       % left=1, right=1
%   [y1, idx] = interp_custom(...)
%
% 输入：
%   x    - 单调递增向量 (N×1)
%   y    - 与 x 等长的函数值向量 (N×1)
%   x1   - 查询点向量
%   left - 左侧外推标志 (0=钳制, 1=线性外推), 默认 0
%   right- 右侧外推标志 (0=钳制, 1=线性外推), 默认 0
%
% 输出：
%   y1  - 与 x1 同尺寸的插值/外推结果
%   idx - 区间索引, 满足 x(idx-1) < x1 <= x(idx)

    % 默认参数
    if nargin < 4, left = 0; end
    if nargin < 5, right = 0; end

    % 查找区间索引
    idx = lookup(x, x1);

    % 两侧端点和对应值
    % xl = x(idx - 1);
    % xr = x(idx);
    % yl = y(idx - 1);
    % yr = y(idx);

    % 线性插值
    % y1 = yl + (yr - yl) ./ (xr - xl) .* (x1(:) - xl);
    y1 = interp1(x, y, x1, 'linear', 'extrap');

    % 左侧处理
    below = x1 < x(1);
    if left == 0
        y1(below) = y(1);
    % else
    %     % 在 [x(1), x(2)] 基础上外推
    %     slopeL = (y(2) - y(1)) / (x(2) - x(1));
    %     y1(below) = y(1) + slopeL .* (x1(below) - x(1));
    end

    % 右侧处理
    above = x1 > x(end);
    if right == 0
        y1(above) = y(end);
    % else
    %     % 在 [x(end-1), x(end)] 基础上外推
    %     slopeR = (y(end) - y(end-1)) / (x(end) - x(end-1));
    %     y1(above) = y(end) + slopeR .* (x1(above) - x(end));
    end
end

function [Va, g, c] = EGMStepBack(Va_p, Xt, par, grid)
% EGMStepBack  内生格点法回推一步
%   给定下一期资产边际价值 Va_p(t+1)（na×ne 矩阵）和当期价格 Xt = [W, R]，
%   输出当期资产边际价值 Va（na×ne），储蓄 g（na×ne），消费 c（na×ne）策略。
%
% 输入：
%   Va_p - (na×ne) 矩阵, grid.a 网格上对每个 shock 的资产边际价值 V'(b')
%   Xt   - 长度 2 向量 [W, R]
%   par  - 参数结构体, 包含 beta, gamma
%   grid - 网格结构体, 包含:
%            a    - 资产网格 (na×1)
%            e    - 收入 shock 向量 (ne×1)
%            Πe   - 转移矩阵 (ne×ne)
%
% 输出：
%   Va - 当期资产边际价值 (na×ne)
%   g  - 当期储蓄决策 b'(a,e) (na×ne)
%   c  - 当期消费决策 c(a,e) (na×ne)

    % 解包参数
    beta  = par.beta;
    gamma = par.gamma;
    W     = Xt(1);
    R     = Xt(2);

    % 1. 当期边际效用 u'(c) = beta * E[Va_p(b',e')]
    uc = beta * Va_p * grid.Pi_e.';  % (na×ne)

    % 2. 逆边际效用求消费: c_implied = (u'^(-1))
    c_implied = uPrimeInv(par, uc);  

    % 3. 预算约束: 资产 a = (b' + c_implied - W*e_j) / R
    % grid.e' 是 1×ne, 扩展到 na×ne
    currentAssets = (grid.a + c_implied - W * grid.e.') / R;  % (na×ne)

    % 4. 插值: 将 b'(a,e) 策略映射回原资产网格
    na = grid.na;
    ne =grid.ne;
    g  = zeros(na, ne);
    for j = 1:ne
        % 对每列 currentAssets(:,j) 在 grid.a 上线性外推
        g(:,j) = interp_custom(currentAssets(:,j), grid.a, grid.a, 0, 1);
        % 这里插值近似的是状态转移方程
    end

    % 5. 借贷约束: b' >= a_min (grid.a(1))
    g = max(g, grid.a(1));

    % 6. 消费 c(a,e) = R*a + W*e - b'
    c = R * grid.a + W * grid.e.' - g;

    % 7. 当期资产边际价值 Va = R * u'(c)
    Va = R * uPrime(par, c);   % uPrime 用户需定义
end

function [Va, g, c] = SolveEGM(c0, Xt, par, grid)
% SOLVEEGM  迭代内生格点法至收敛
% 输入:
%   c0   - 初始消费矩阵 (na×ne)
%   Xt   - 当期价格向量 [W, R]
%   par  - 参数结构体，包含 beta, gamma
%   grid - 网格结构体，包含 a, e, Πe
% 输出:
%   Va - 收敛后的资产边际价值 (na×ne)
%   g  - 储蓄策略矩阵 (na×ne)
%   c  - 消费策略矩阵 (na×ne)

    Va   = Xt(2)*uPrime(par, c0);
    tol   = 1e-10;
    maxit = 10000;

    for it = 1:maxit
        [Va1, g1, c1] = EGMStepBack(Va, Xt, par, grid);

        if mod(it,50) == 0
            test = abs(Va1 - Va) ./ (abs(Va) + tol);
            if all(test(:) < tol)
                Va = Va1;
                break;
            end
        end
        Va = Va1;

    end
    [Va, g, c] = EGMStepBack(Va, Xt, par, grid);
end

