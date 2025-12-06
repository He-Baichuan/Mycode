clear;
clc;


param = SettingPar_SSJ;
mY = repmat(param.y,1,param.N_a);
r = 0.01/4;
[Va, ap, cp] = Policy_SS(mY, r, param);
% 假设已有变量：
% c - 矩阵，形状为 [n_y, n_a]，其中 n_y 是收入水平数，n_a 是资产网格点数
% a_grid - 行向量，长度为 n_a，表示资产网格
% r - 标量，利率
% a - 矩阵，形状与 c 相同，表示每个点对应的资产值

% 初始化 mpcs 矩阵
[n_y, n_a] = size(cp);
mpcs = zeros(n_y, n_a);  % 等价于 np.empty_like(c)
   % for i = 1:n_y
   %      % gradient计算梯度，第二个参数是间距
   %      dcdx = gradient(cp(i, :), param.gridA);
   %      mpcs(i, :) = dcdx / (1 + r);
   %  end
% 
% 1. 内部点：对称差分
% Python: mpcs[:, 1:-1] = (c[:, 2:] - c[:, 0:-2]) / (a_grid[2:] - a_grid[:-2]) / (1+r)
mpcs(:, 2:end-1) = (cp(:, 3:end) - cp(:, 1:end-2)) ./ ...
                  ((param.gridA(3:end,1) - param.gridA(1:end-2,1)) ./ (1 + r))';

% 2. 边界点：非对称差分
% 左边界
% Python: mpcs[:, 0] = (c[:, 1] - c[:, 0]) / (a_grid[1] - a_grid[0]) / (1+r)
mpcs(:, 1) = (cp(:, 2) - cp(:, 1)) / ((param.gridA(2) - param.gridA(1)) / (1 + r))';

% 右边界
% Python: mpcs[:, -1] = (c[:, -1] - c[:, -2]) / (a_grid[-1] - a_grid[-2]) / (1+r)
mpcs(:, end) = (cp(:, end) - cp(:, end-1)) / ((param.gridA(end) - param.gridA(end-1)) / (1 + r))';

% 3. 特殊情况：当资产等于网格最小值时，MPC=1
% Python: mpcs[a == a_grid[0]] = 1
mpcs(ap == param.gridA(1)) = 1;

figure
plot(param.gridA,cp)
ylim([0 1])
xlim([0 1])
% legend('y1','y2','y3','y4','y5','y6','y7')
legend_label1 = arrayfun(@(ye) sprintf('y=%.2f', ye), param.y, 'UniformOutput', false);
legend(legend_label1);

figure
plot(param.gridA,mpcs)
ylim([0 1])
xlim([0 0.3])
% legend('y1','y2','y3','y4','y5','y6','y7')
% legend_label2 = arrayfun(@(i) sprintf('MPC Curve %d', i), 1:size(mpcs,1), 'UniformOutput', false);
legend(legend_label1);


function param = SettingPar_SSJ()
rho = 0.975;
sigma = 0.7;
N_e = 7;
beta = 1-0.08/4;
gamma = 1;
a_min = 0;
a_max = 10000;
nu = 0.025;
N_a = 500;
if nu == 0.0
        % 对应 Julia: collect(range(-ϕ, a_max, nA))
        gridA = linspace(a_min, a_max, N_a).';
else
        % 对应 Julia: step = 0:(nA-1); 非线性拉伸网格
        step  = 0:(N_a-1);
        gridA = a_min + (a_max - a_min) .* (( (1 + nu).^step - 1.0) ./ ((1 + nu)^(N_a - 1) - 1.0));
        gridA = gridA(:);   % 转成列向量
end
[~, transZ] = rouwenhorst(N_e,rho, sigma );
inv_distZ = 1/N_e *ones(N_e,1);
for t = 1:10000
    inv_distZ_new = inv_distZ' * transZ;
    crit = max(abs(inv_distZ_new - inv_distZ'));
    inv_distZ = inv_distZ_new';
    if crit<1e-10
        break
    end
end

gridZ = (0:N_e - 1)';
alpha  = 2*sigma/sqrt(N_e-1);
gridZ = alpha*gridZ;
gridZ = exp(gridZ);
y = gridZ/(sum(gridZ.*inv_distZ));
param.y = y;
param.beta = beta;
param.gamma = gamma;
param.gridA = gridA;
param.transZ = transZ;
param.inv_distZ = inv_distZ;
param.N_e = N_e;
param.N_a = N_a;
end

function [Va_new, ap_new, cp] =  StepBackEGM(Va, y, r, param )
gridA = param.gridA;
transZ = param.transZ;
beta = param.beta;
gamma = param.gamma;
N_e = param.N_e;
N_a = param.N_a;
gAA = repmat(gridA',N_e,1);
RHS = beta * transZ * Va;
c = RHS.^(-gamma);
A_endo = (c + gAA - y)/(1+r);
ap_new = zeros(N_e, N_a);
coh = y+(1+r)*gAA;
        for i = 1:N_e
            % A_endo(:,j) : 内生格点 (当前资产 a 对应的“可支配收入”)
            % gridA       : 在 Julia 里是 ygrid，这里 y = a' 的值
            % gridA       : xtarget, 目标是要在外生资产网格上得到 a'(a,s)

            % 基本的线性插值（允许外推，后面自己截断）：
            ap_temp = interp1( A_endo(i,:), ...   % xgrid
                            gridA', ...          % ygrid
                            gridA', ...          % xquery
                            'linear', 'extrap');

            % 处理借款约束：小于下界的一律截到 gridA(1)
            ap_temp(ap_temp < gridA(1)) = gridA(1);

            % 也可以视情况对上界做截断：
            % ap_temp(ap_temp > gridA(end)) = gridA(end);

            ap_new(i,:) = ap_temp(:)';
        end
cp = coh - ap_new;
Va_new = (1+r)*cp.^(-1/gamma);
end

function [Va, ap, cp] = Policy_SS(y, r, param)
gridA = param.gridA;
N_e = param.N_e;
gamma = param.gamma;

mAsset = repmat(gridA',N_e,1);
coh = y + (1 + r)*mAsset;
c = coh*0.5;
Va = c.^(-1/gamma) * (1+r);
    
    for t = 1:10000
        [Va_new, ap, cp] = StepBackEGM(Va,y,r,param);
        if max(abs(Va_new - Va))<1e-8
            break
        end
        Va = Va_new;
    end

end
