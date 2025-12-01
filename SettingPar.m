function par = SettingPar()
%  SettingPar
%  生成转移动态所需的参数与各类序列
%  （这是 Julia 版本 SettingPar_TransitionDynamics 的直接翻译）

    %% ----------- 标准参数（对应 Julia 中的默认 keyword 参数） -----------
    nA   = 300;
    nS   = 7;
    alpha = 0.33;
    beta  = 0.96;
    gamma = 2.0;
    delta = 0.05;
    phi   = 0.0;
    rho   = 0.9;
    sigma = 0.03;
    Zbar  = 0.0;   % 总冲击的 log
    rho_z = 0.95;
    T     = 150;

    % 税收与政府支出相关参数
    ShockSize = 0.01;
    tau_l0 = 0.0;   % 第 0 期劳动税率
    tau_l1 = 0.0;   % 第 1 期劳动税率
    tau_k0 = 0.0;   % 第 0 期资本税率
    tau_k1 = 0.0;   % 第 1 期资本税率
    G0     = 0.0;   % 第 0 期政府支出比率
    G1     = 0.0;   % 第 1 期政府支出比率

    a_max  = 100.0;
    nu     = 0.025;

    %% ----------- 生成劳动格点（Rouwenhorst） -----------
    % 假定 rouwenhorst 实现为 package 函数：+MyNumerical_method/rouwenhorst.m
    [gridS, transS] = rouwenhorst(nS, rho, sigma);
    gridS = exp(gridS);   % 对数 -> 水平值

    %% ----------- 计算不变分布 inv_S -----------
    inv_S = ones(nS,1) / nS;   % 初始化为均匀分布（列向量）
    for i = 1:10000
        % 对应 Julia: inv_S_new = (inv_S' * transS)'
        inv_S_new = (inv_S' * transS).';
        if max(abs(inv_S_new - inv_S)) < 1e-8
            break;
        end
        inv_S = inv_S_new;
    end
    % 劳动禀赋的稳态均值 Lbar = Σ_s gridS(s) * inv_S(s)
    Lbar = sum(gridS(:) .* inv_S(:));

    %% ----------- 生成资产格点 gridA -----------
    if nu == 0.0
        % 对应 Julia: collect(range(-ϕ, a_max, nA))
        gridA = linspace(-phi, a_max, nA).';
    else
        % 对应 Julia: step = 0:(nA-1); 非线性拉伸网格
        step  = 0:(nA-1);
        gridA = -phi + (a_max + phi) .* (( (1 + nu).^step - 1.0) ./ ((1 + nu)^(nA - 1) - 1.0));
        gridA = gridA(:);   % 转成列向量
    end

    %% ----------- 构造总冲击序列 Shock_Z -----------
    Shock_Z = zeros(T,1);
    Zbar     = exp(Zbar);   % Julia 中先把 log(Zbar) 还原为水平 Zbar

    Shock_Z(1) = log(Zbar);
    Shock_Z(T) = log(Zbar);
    Shock_Z(2) = ShockSize;
    for t = 3:(T-1)
        Shock_Z(t) = rho_z * Shock_Z(t-1);
    end
    Shock_Z = exp(Shock_Z);   % 回到水平值

    %% ----------- 构造政策冲击序列（永久性变化） -----------
    tau_l_seq = tau_l1 * ones(T,1);
    tau_l_seq(1) = tau_l0;

    tau_k_seq = tau_k1 * ones(T,1);
    tau_k_seq(1) = tau_k0;

    Gseq = G1 * ones(T,1);
    Gseq(1) = G0;

    %% ----------- 打包为结构体（对应 Julia 的 named tuple） -----------
    par.alpha   = alpha;
    par.beta    = beta;
    par.gamma   = gamma;
    par.delta   = delta;
    par.phi     = phi;

    par.Zbar    = Zbar;
    par.rho_z   = rho_z;

    par.gridA   = gridA;
    par.gridS   = gridS;
    par.transS  = transS;
    par.Lbar    = Lbar;

    par.T       = T;
    par.nA      = nA;
    par.nS      = nS;
    par.nMeasure = 3;
    par.nAQuadrature = 8;
    par.nSQuadrature = par.nAQuadrature*par.nS;
    par.nMeasureCoefficients = par.nMeasure * par.nS;
    par.Shock_Z    = Shock_Z;
    par.tau_l_seq  = tau_l_seq;
    par.tau_k_seq  = tau_k_seq;
    par.Gseq       = Gseq;

    %（可选额外输出：如果你想在其它地方用到这些基础参数，也可以一并放进去）
    par.rho        = rho;
    par.sigma      = sigma;
    par.a_max      = a_max;
    par.nu         = nu;
    par.ShockSize  = ShockSize;
    par.gridAZeros = scaleDown(par.gridA,-par.phi, par.a_max);
    [vAssetsGridQuadratureZeros,vQuadratureWeights] = computeGaussLegendreQuadrature(par.nAQuadrature);
    par.gridAQuadrature = scaleUp(vAssetsGridQuadratureZeros,-par.phi + 0.1,par.a_max);
    par.QuadratureWeights = vQuadratureWeights;
    par.inv_S = inv_S;
end
