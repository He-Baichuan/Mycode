clear; clc;

% 主程序顺序：
% 1. 初始化 Fortran 模型参数；
% 2. 用 CEtools 建立资产维度的 collocation basis；
% 3. 按年龄倒推求解 Mongey Case I 的 (V, EV) 系数；
% 4. 用 Mongey 的 QX/Q 矩阵推进生命周期分布；
% 5. 聚合并输出诊断结果。

tic
par = initialize_model();
col = initialize_collocation(par);
sol = solve_collocation(par, col);
dist = get_distribution_mongey(sol, par, col);
agg = aggregation_mongey(sol, dist, par, col);
toc
output_results(sol, dist, agg, par, col);

disp('Collocation household solver and Mongey distribution finished.');

function par = initialize_model()
    % 参数基本沿用 Fortran 的 globals 和 initialize 子程序。
    par.JJ = 80; % 年龄
    par.JR = 45; % 退休年龄
    par.NA = 200; % 资产网格数
    par.NP = 2; % 生产率状态数
    par.NS = 7; % 持久性冲击状态数

    par.gamma = 0.50;
    par.egam = 1 - 1/par.gamma;
    par.beta = 0.98;
    par.nu = 0.335;

    par.sigma_theta = 0.242;
    par.sigma_eps = 0.022;
    par.rho = 0.985;

    par.a_l = 0.0;
    par.a_u = 200.0;
    par.a_grow = 0.05;
    par.c_min = 1e-10;
    par.l_min = 0.0;
    par.l_max = 1.0 - 1e-10;

    par.r = 0.04; % 外生利率
    par.w = 1.0; % 外生工资率

    par.psi = [ ...
        1.00000, 0.99923, 0.99914, 0.99914, 0.99912, ...
        0.99906, 0.99908, 0.99906, 0.99907, 0.99901, ...
        0.99899, 0.99896, 0.99893, 0.99890, 0.99887, ...
        0.99886, 0.99878, 0.99871, 0.99862, 0.99853, ...
        0.99841, 0.99835, 0.99819, 0.99801, 0.99785, ...
        0.99757, 0.99735, 0.99701, 0.99676, 0.99650, ...
        0.99614, 0.99581, 0.99555, 0.99503, 0.99471, ...
        0.99435, 0.99393, 0.99343, 0.99294, 0.99237, ...
        0.99190, 0.99137, 0.99085, 0.99000, 0.98871, ...
        0.98871, 0.98721, 0.98612, 0.98462, 0.98376, ...
        0.98226, 0.98062, 0.97908, 0.97682, 0.97514, ...
        0.97250, 0.96925, 0.96710, 0.96330, 0.95965, ...
        0.95619, 0.95115, 0.94677, 0.93987, 0.93445, ...
        0.92717, 0.91872, 0.91006, 0.90036, 0.88744, ...
        0.87539, 0.85936, 0.84996, 0.82889, 0.81469, ...
        0.79705, 0.78081, 0.76174, 0.74195, 0.72155, ...
        0.00000]';

    par.eff = zeros(par.JJ, 1);
    par.eff(1:par.JR-1) = [ ...
        1.0000, 1.0719, 1.1438, 1.2158, 1.2842, 1.3527, ...
        1.4212, 1.4897, 1.5582, 1.6267, 1.6952, 1.7217, ...
        1.7438, 1.7748, 1.8014, 1.8279, 1.8545, 1.8810, ...
        1.9075, 1.9341, 1.9606, 1.9623, 1.9640, 1.9658, ...
        1.9675, 1.9692, 1.9709, 1.9726, 1.9743, 1.9760, ...
        1.9777, 1.9700, 1.9623, 1.9546, 1.9469, 1.9392, ...
        1.9315, 1.9238, 1.9161, 1.9084, 1.9007, 1.8354, ...
        1.7701, 1.7048]';

    par.pen = zeros(par.JJ, 1);
    par.pen(par.JR:par.JJ) = 0.5*sum(par.eff)/(par.JR-1)*0.33;

    par.dist_theta = ones(par.NP, 1)/par.NP;
    par.theta = [exp(-sqrt(par.sigma_theta)); exp(sqrt(par.sigma_theta))];
    % Fortran toolbox 的 discretize_AR 使用 Rouwenhorst，不是 Tauchen。
    [par.eta, par.pi] = rouwenhorst_local(par.rho, 0.0, par.sigma_eps, par.NS);
    par.eta = exp(par.eta(:));

    % 资产网格完全复刻 toolbox.f90 的 grid_Cons_Grow。
    par.a = grid_cons_grow(par.a_l, par.a_u, par.a_grow, par.NA);
    par.a_bor = zeros(par.JJ, par.NP);

    % 预先计算工资率和非劳动现金在手，避免 Newton/黄金分割内部重复计算。
    par.wage = zeros(par.JJ, par.NP, par.NS);
    for j = 1:par.JJ
        for ip = 1:par.NP
            for is = 1:par.NS
                if j < par.JR
                    par.wage(j, ip, is) = ...
                        par.w*par.eff(j)*par.theta(ip)*par.eta(is);
                end
            end
        end
    end

    par.cash_on_hand = zeros(par.JJ, par.NA+1, par.NP, par.NS);
    for j = 1:par.JJ
        for ia = 1:par.NA+1
            for ip = 1:par.NP
                for is = 1:par.NS
                    par.cash_on_hand(j, ia, ip, is) = ...
                        (1+par.r)*par.a(ia) + par.pen(j);
                end
            end
        end
    end
end

function col = initialize_collocation(par)
    % 资产维度使用 CEtools 的线性 spline basis。breakpoints 直接取 Fortran 资产网格，
    % 因此 collocation nodes 与 par.a 一一对应。
    col.fspace_a = fundef({'spli', par.a, 0, 1});
    col.a_nodes = funnode(col.fspace_a);
    col.Phi_a = funbas(col.fspace_a, col.a_nodes);
    col.nx = numel(col.a_nodes); % 资产格点数
    col.ns = par.NS;
    col.m = col.nx*col.ns;  % 状态变量组合数

    % Mongey 堆叠顺序：资产快变、冲击慢变，故 E[V] = kron(pi,I_a)*V。
    col.Phi = kron(speye(par.NS), sparse(col.Phi_a)); % 配点法下价值函数对应的基 (basis)
    col.PEV = kron(sparse(par.pi), speye(col.nx)); % 配点法下期望价值函数对应的基 (basis)

    % 分布网格先取求解网格。若以后需要更细分布网格，只需替换这里。
    col.fspace_d = col.fspace_a;
    col.a_dist = col.a_nodes;
    col.nxd = col.nx;
    col.nd = col.nxd*par.NS;

    if col.nx ~= par.NA + 1
        error('Unexpected number of collocation nodes: %d', col.nx);
    end
end

function sol = solve_collocation(par, col)
    sol = repmat(empty_solution(par, col), par.JJ, 1);

    % 终期没有 continuation，全部资源用于消费，aplus=0。
    j = par.JJ; % 年龄 ---> 生产率 ---> 持久性冲击 
    for ip = 1:par.NP
        Vmat = zeros(col.nx, par.NS);
        Cmat = zeros(col.nx, par.NS);
        Lmat = zeros(col.nx, par.NS);
        for is = 1:par.NS
            resources = squeeze(par.cash_on_hand(j, :, ip, is));
            resources = resources(:);
            Cmat(:, is) = max(resources, par.c_min);
            Vmat(:, is) = utility(Cmat(:, is), Lmat(:, is), par); % 先计算最后一期的价值函数
        end
        EVmat = reshape(col.PEV*Vmat(:), col.nx, par.NS);
        sol(j).coefV(:, ip, :) = reshape(col.Phi\Vmat(:), col.nx, 1, par.NS);
        sol(j).coefEV(:, ip, :) = reshape(col.Phi\EVmat(:), col.nx, 1, par.NS);
        sol(j).V(:, ip, :) = reshape(Vmat, col.nx, 1, par.NS);
        sol(j).EV(:, ip, :) = reshape(EVmat, col.nx, 1, par.NS);
        sol(j).c(:, ip, :) = reshape(Cmat, col.nx, 1, par.NS);
        sol(j).l(:, ip, :) = reshape(Lmat, col.nx, 1, par.NS);
        sol(j).aplus(:, ip, :) = 0;
    end
    sol(j).residual = 0;
    sol(j).flag = 0;
    sol(j).iterations = 0;

    for j = par.JJ-1:-1:1
        fprintf('Solving age %3d\n', j);
        for ip = 1:par.NP
            % 固定效应 ip 无转移，因此每个 ip 单独解一个 NX*NS 系统。
            % 给定下一期 EV 系数，先在各节点用 goldenx 求 Bellman RHS 和 policy。
            [rhsV, c_policy, l_policy, a_policy] = policy_rhs_age_ip(j, ip, sol(j+1), par, col);

            z0 = initial_guess(sol(j+1), ip, par, col);
            ctx.Phi = col.Phi;
            ctx.PEV = col.PEV;
            ctx.rhsV = rhsV;
            ctx.jac = [col.Phi, sparse(col.m, col.m); ...
                       -col.PEV*col.Phi, col.Phi];

            [z, fval, flag, it] = newton(@obj_newton, z0, ctx);
            [coefV, coefEV] = unpack_coefficients(z, col);

            Vvec = col.Phi*coefV(:);
            EVvec = col.Phi*coefEV(:);

            sol(j).coefV(:, ip, :) = reshape(coefV, col.nx, 1, par.NS);
            sol(j).coefEV(:, ip, :) = reshape(coefEV, col.nx, 1, par.NS);
            sol(j).V(:, ip, :) = reshape(reshape(Vvec, col.nx, par.NS), col.nx, 1, par.NS);
            sol(j).EV(:, ip, :) = reshape(reshape(EVvec, col.nx, par.NS), col.nx, 1, par.NS);
            sol(j).c(:, ip, :) = reshape(c_policy, col.nx, 1, par.NS);
            sol(j).l(:, ip, :) = reshape(l_policy, col.nx, 1, par.NS);
            sol(j).aplus(:, ip, :) = reshape(a_policy, col.nx, 1, par.NS);
            sol(j).residual = max(sol(j).residual, norm(fval, inf));
            sol(j).flag = max(sol(j).flag, flag);
            sol(j).iterations = max(sol(j).iterations, it);
        end
    end
end

function s = empty_solution(par, col)
    % 每个年龄下的系数
    s.coefV = zeros(col.nx, par.NP, par.NS);
    s.coefEV = zeros(col.nx, par.NP, par.NS);
    s.V = zeros(col.nx, par.NP, par.NS);
    s.EV = zeros(col.nx, par.NP, par.NS);
    s.c = zeros(col.nx, par.NP, par.NS);
    s.l = zeros(col.nx, par.NP, par.NS);
    s.aplus = zeros(col.nx, par.NP, par.NS);
    s.residual = 0;
    s.flag = 0;
    s.iterations = 0;
end

function z0 = initial_guess(sol_next, ip, par, col)
    coefV0 = reshape(sol_next.coefV(:, ip, :), col.nx, par.NS); % 维度 nA x nS
    coefEV0 = reshape(sol_next.coefEV(:, ip, :), col.nx, par.NS); % 维度 nA x nS
    z0 = [coefV0(:); coefEV0(:)];
end

function [res, jac] = obj_newton(z, ctx)
    % Mongey Case I 残差：
    % g1 = Phi*c  - [u(c_policy)+beta*psi*EV_next(aplus)]
    % g2 = Phi*ce - (pi kron I_a)*Phi*c
    % 其中 ce 是 expected continuation value 的系数，不是 Fortran 的 certainty equivalent。
    m = numel(ctx.rhsV);
    c = z(1:m);
    ce = z(m+1:end);

    g1 = ctx.Phi*c - ctx.rhsV;
    g2 = ctx.Phi*ce - ctx.PEV*ctx.Phi*c;
    res = [g1; g2];

    if nargout > 1
        jac = ctx.jac;
    end
end

function [rhsV, c_policy, l_policy, a_policy] = policy_rhs_age_ip(j, ip, sol_next, par, col)
    % 对当前年龄 j 和固定效应 ip，在所有资产节点和 eta 状态上求单点储蓄政策。
    % goldenx 是向量化黄金分割，可一次处理同一个 is 下的全部资产节点。
    rhsV = zeros(col.nx, par.NS);
    c_policy = zeros(col.nx, par.NS);
    l_policy = zeros(col.nx, par.NS);
    a_policy = zeros(col.nx, par.NS);

    for is = 1:par.NS
        resources = squeeze(par.cash_on_hand(j, :, ip, is));
        resources = resources(:);
        wage = par.wage(j, ip, is);
        ap_low = zeros(col.nx, 1);
        ap_high = resources + wage*par.l_max - par.c_min;
        ap_high = min(par.a_u, max(0, ap_high));

        coefEV_next = reshape(sol_next.coefEV(:, ip, is), col.nx, 1);
        obj = @(ap) bellman_rhs(ap, resources, coefEV_next, j, ip, is, par, col);
        [ap, val] = goldenx(obj, ap_low, ap_high);

        [cons, lab] = choices_from_aplus(ap, resources, j, ip, is, par);
        a_policy(:, is) = ap;
        c_policy(:, is) = cons;
        l_policy(:, is) = lab;
        rhsV(:, is) = val;
    end
    rhsV = rhsV(:);
end

function val = bellman_rhs(aplus, resources, coefEV_next, j, ip, is, par, col) 
    % Bellman 右端使用 Mongey 的普通 expected value：
    % u(c) + beta*psi(j+1)*EV_{j+1}(aplus,ip,is)。
    [cons, lab] = choices_from_aplus(aplus, resources, j, ip, is, par);
    cont = funeval(coefEV_next, col.fspace_a, aplus); % 期望价值函数 EV(a+)
    cont = reshape(cont, size(aplus));
    val = utility(cons, lab, par) + par.beta*par.psi(j+1)*cont;
end

function [cons, lab] = choices_from_aplus(aplus, resources, j, ip, is, par)
    wage = par.wage(j, ip, is);
    if j < par.JR && wage > 0
        lab = par.nu + (1-par.nu)*(aplus - resources)/wage;
        lab = min(max(lab, par.l_min), par.l_max);
    else
        lab = zeros(size(aplus));
    end
    cons = max(resources + wage*lab - aplus, par.c_min);
end

function [coefV, coefEV] = unpack_coefficients(z, col)
    m = col.m;
    coefV = reshape(z(1:m), col.nx, col.ns);
    coefEV = reshape(z(m+1:end), col.nx, col.ns);
end

function u = utility(c, lab, par)
    c = max(c, par.c_min);
    lab = min(max(lab, par.l_min), par.l_max);
    flow = c.^par.nu .* (1-lab).^(1-par.nu);
    if abs(par.egam) < 1e-12
        u = log(flow);
    else
        u = flow.^par.egam/par.egam;
    end
end

function a = grid_cons_grow(left, right, growth, NA)
    % CE-Fortran toolbox.f90: grid_Cons_Grow。
    n = NA;
    h = (right-left)/((1+growth)^n - 1);
    a = h*((1+growth).^(0:n)' - 1) + left;
end

function [z, pi] = rouwenhorst_local(rho, mu, sigma_eps, NS)
    % 复刻 toolbox.f90 的 discretize_AR：
    % sigma_eps 在原注释里称 shock variance，但代码按该数值直接计算 sigma_eta。
    sigma_eta = sigma_eps/(1-rho^2);
    psi = sqrt(NS-1)*sqrt(sigma_eta);
    z = linspace(-psi, psi, NS)' + mu;
    pi = rouwenhorst_matrix_local(rho, NS);
end

function pi = rouwenhorst_matrix_local(rho, n)
    % 递归 Rouwenhorst 转移矩阵，内部行除以 2 与 Fortran 保持一致。
    p = (1+rho)/2;
    if n == 2
        pi = [p, 1-p; 1-p, p];
        return
    end

    pi_old = rouwenhorst_matrix_local(rho, n-1);
    pi = zeros(n, n);
    pi(1:n-1, 1:n-1) = pi(1:n-1, 1:n-1) + p*pi_old;
    pi(1:n-1, 2:n) = pi(1:n-1, 2:n) + (1-p)*pi_old;
    pi(2:n, 1:n-1) = pi(2:n, 1:n-1) + (1-p)*pi_old;
    pi(2:n, 2:n) = pi(2:n, 2:n) + p*pi_old;
    pi(2:n-1, :) = pi(2:n-1, :)/2;
end

function dist = get_distribution_mongey(sol, par, col)
    % Mongey 分布对象 L(s_d)：行状态按 (a_1,is_1),...,(a_N,is_1),(a_1,is_2),...
    % 生命周期模型中 Q_j 随年龄变化，所以逐年龄推进 L_{j+1}=Q_j'*L_j，
    % 不求 stationary eigenvector。
    dist.L = zeros(col.nd, par.NP, par.JJ);
    dist.mass = zeros(par.JJ, par.NP);
    dist.mass_error = zeros(par.JJ, 1);
    dist.q_rowsum_error = zeros(par.JJ-1, par.NP);
    dist.is_initial = ceil(par.NS/2);

    for ip = 1:par.NP
        % 初始 cohort 位于 eta 中间状态，初始资产为 Fortran 的 a=0。
        % 用分布网格的线性 spline basis 把该点质量分配到相邻资产点。
        a0 = col.a_dist(1);
        qx0 = funbas(col.fspace_d, a0); % young 的 lottery 方法
        rows0 = (dist.is_initial-1)*col.nxd + (1:col.nxd);
        dist.L(rows0, ip, 1) = qx0(:)*par.dist_theta(ip); % 生产率初始（不变）分布 ---> 资产的分配权重
    end

    dist.mass(1, :) = squeeze(sum(dist.L(:, :, 1), 1)); % mass_{j} = sum(L) = dist(eta) * sum(q) = dist(eta)
    dist.mass_error(1) = abs(sum(dist.mass(1, :)) - 1); % sum(mass_{j}) = sum(dist(eta)) = 1

    for j = 1:par.JJ-1
        for ip = 1:par.NP
            % QX 是 Mongey notes 中的资产转移矩阵：每一行是 aplus(s)
            % 在分布资产网格上的线性插值权重。
            ap = reshape(sol(j).aplus(:, ip, :), col.nxd, par.NS);
            ap = min(max(ap(:), col.a_dist(1)), col.a_dist(end)); % 限定在支撑集里面
            QX = sparse(funbas(col.fspace_d, ap));

            % Q 同时包含内生资产转移和外生 eta 转移。固定效应 ip 不转移。
            Q = build_transition_mongey(QX, par.pi, col.nxd, par.NS);

            dist.q_rowsum_error(j, ip) = full(max(abs(sum(Q, 2) - 1))); % 验证行和为1
            dist.L(:, ip, j+1) = Q' * dist.L(:, ip, j);
        end
        dist.mass(j+1, :) = squeeze(sum(dist.L(:, :, j+1), 1));
        dist.mass_error(j+1) = abs(sum(dist.mass(j+1, :)) - 1);
    end
end

function Q = build_transition_mongey(QX, pi, nxd, NS)
    % 显式实现 Mongey 的 Q=dprod(QX,QZ)。
    % 对当前状态 (a,is)，下一期列块 isp 的资产权重为 pi(is,isp)*QX(row,:)。
    % 持久性冲击的演化对应的是年龄的增加
    nd = nxd*NS;
    Q = spalloc(nd, nd, nnz(QX)*NS);
    for is = 1:NS
        rows = (is-1)*nxd + (1:nxd);
        for isp = 1:NS
            cols = (isp-1)*nxd + (1:nxd);
            Q(rows, cols) = pi(is, isp)*QX(rows, :);
        end
    end
end

function agg = aggregation_mongey(sol, dist, par, col)
    % 所有 cohort 统计量都直接用 L(s_d)'*policy_or_state_vector 计算。
    % 这里的资产为 Fortran 口径 a。
    agg.c_coh = zeros(par.JJ, 1);
    agg.y_coh = zeros(par.JJ, 1);
    agg.l_coh = zeros(par.JJ, 1);
    agg.h_coh = zeros(par.JJ, 1);
    agg.income_coh = zeros(par.JJ, 1);
    agg.pen_coh = zeros(par.JJ, 1);
    agg.a_coh = zeros(par.JJ, 1);
    agg.v_coh = zeros(par.JJ, 1);
    agg.cv_c = zeros(par.JJ, 1);
    agg.cv_y = zeros(par.JJ, 1);
    agg.cv_l = zeros(par.JJ, 1);
    agg.cv_h = zeros(par.JJ, 1);
    agg.corr_hl = zeros(par.JJ, 1);
    agg.frac_bor = zeros(par.JJ, 1);

    for j = 1:par.JJ
        c2 = 0;
        y2 = 0;
        l2 = 0;
        h2 = 0;
        hl = 0;
        for ip = 1:par.NP
            L = dist.L(:, ip, j);
            % reshape 后再 (:) 可保持 Mongey 的资产快变、冲击慢变堆叠顺序。
            cvec = reshape(sol(j).c(:, ip, :), col.nxd, par.NS);
            lvec = reshape(sol(j).l(:, ip, :), col.nxd, par.NS);
            avec = repmat(col.a_dist, par.NS, 1);
            vvec = reshape(sol(j).V(:, ip, :), col.nxd, par.NS);
            apvec = reshape(sol(j).aplus(:, ip, :), col.nxd, par.NS);

            hmat = zeros(col.nxd, par.NS);
            ymat = zeros(col.nxd, par.NS);
            for is = 1:par.NS
                hmat(:, is) = par.eff(j)*par.theta(ip)*par.eta(is);
                ymat(:, is) = hmat(:, is).*lvec(:, is);
            end

            cvec = cvec(:);
            lvec = lvec(:);
            hvec = hmat(:);
            yvec = ymat(:);
            vvec = vvec(:);
            apvec = apvec(:);

            agg.c_coh(j) = agg.c_coh(j) + cvec' * L;
            agg.y_coh(j) = agg.y_coh(j) + yvec' * L;
            agg.l_coh(j) = agg.l_coh(j) + lvec' * L;
            agg.h_coh(j) = agg.h_coh(j) + hvec' * L;
            agg.a_coh(j) = agg.a_coh(j) + avec' * L;
            agg.v_coh(j) = agg.v_coh(j) + vvec' * L;
            agg.frac_bor(j) = agg.frac_bor(j) + double(apvec < 1e-6)' * L;

            c2 = c2 + (cvec.^2)' * L;
            y2 = y2 + (yvec.^2)' * L;
            l2 = l2 + (lvec.^2)' * L;
            h2 = h2 + (hvec.^2)' * L;
            hl = hl + (hvec.*lvec)' * L;
        end
        agg.pen_coh(j) = par.pen(j)*sum(dist.L(:, :, j), 'all');
        agg.income_coh(j) = agg.y_coh(j) + par.r*agg.a_coh(j);
        agg.cv_c(j) = sqrt(max(c2 - agg.c_coh(j)^2, 0))/max(agg.c_coh(j), par.c_min);
        agg.cv_y(j) = sqrt(max(y2 - agg.y_coh(j)^2, 0))/max(agg.y_coh(j), par.c_min);
        agg.cv_l(j) = sqrt(max(l2 - agg.l_coh(j)^2, 0))/max(agg.l_coh(j), par.c_min);
        agg.cv_h(j) = sqrt(max(h2 - agg.h_coh(j)^2, 0))/max(agg.h_coh(j), par.c_min);
        denom = sqrt(max(h2 - agg.h_coh(j)^2, 0))*sqrt(max(l2 - agg.l_coh(j)^2, 0));
        agg.corr_hl(j) = (hl - agg.h_coh(j)*agg.l_coh(j))/max(denom, par.c_min);
    end
    agg.frac_bor(par.JJ) = 1;
end

function output_results(sol, dist, agg, par, col)
    % 输出两类文件：policy/Newton 诊断，以及包含分布推进后的 cohort aggregation。
    ages = 20 + (1:par.JJ)';
    residual = zeros(par.JJ, 1);
    min_c = zeros(par.JJ, 1);
    min_l = zeros(par.JJ, 1);
    max_l = zeros(par.JJ, 1);
    max_aplus = zeros(par.JJ, 1);
    upper_share = zeros(par.JJ, 1);
    ev_error = zeros(par.JJ, 1);

    for j = 1:par.JJ
        residual(j) = sol(j).residual;
        min_c(j) = min(sol(j).c(:));
        min_l(j) = min(sol(j).l(:));
        max_l(j) = max(sol(j).l(:));
        max_aplus(j) = max(sol(j).aplus(:));
        upper_share(j) = mean(sol(j).aplus(:) > par.a_u - 1e-6);
        ev_error(j) = max_ev_error(sol(j), par, col);
    end

    fid = fopen('output_collocation_labor_policy.out', 'w');
    if fid < 0
        error('Could not open output_collocation_labor_policy.out for writing.');
    end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, ' AGE    RESIDUAL       MIN_C       MIN_L       MAX_L      MAX_APLUS   UPPER_SHARE   EV_ERROR\n');
    for j = 1:par.JJ
        fprintf(fid, '%4d %12.4e %12.4e %11.4e %11.4e %12.4f %12.4e %12.4e\n', ...
            ages(j), residual(j), min_c(j), min_l(j), max_l(j), max_aplus(j), upper_share(j), ev_error(j));
    end
    clear cleanup

    fid = fopen('output_collocation_labor_distribution.out', 'w');
    if fid < 0
        error('Could not open output_collocation_labor_distribution.out for writing.');
    end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, ' AGE       MASS_ERR      CONS     HOURS   EARNINGS     INCOME       PENS     ASSETS      CV_C      CV_L      CV_Y      CV_H   CORR_HL  FRAC_BOR\n');
    for j = 1:par.JJ
        fprintf(fid, '%4d %12.4e %9.3f %9.3f %10.3f %10.3f %10.3f %10.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n', ...
            ages(j), dist.mass_error(j), agg.c_coh(j), agg.l_coh(j), agg.y_coh(j), agg.income_coh(j), ...
            agg.pen_coh(j), agg.a_coh(j), agg.cv_c(j), agg.cv_l(j), agg.cv_y(j), agg.cv_h(j), agg.corr_hl(j), agg.frac_bor(j));
    end
    clear cleanup

    save('prog10_02_collocation_solution.mat', 'sol', 'dist', 'agg', 'par', 'col');
    make_policy_plots(sol, par, col, ages, residual);
    make_distribution_plots(agg, dist, ages);
end

function err = max_ev_error(sj, par, col)
    % 检查 Mongey 第二个残差方程 EV=(pi kron I_a)*V。
    err = 0;
    for ip = 1:par.NP
        Vmat = reshape(sj.V(:, ip, :), col.nx, par.NS);
        EVmat = reshape(sj.EV(:, ip, :), col.nx, par.NS);
        target = reshape(col.PEV*Vmat(:), col.nx, par.NS);
        err = max(err, max(abs(EVmat(:) - target(:))));
    end
end

function make_policy_plots(sol, par, col, ages, residual)
    % 保存若干 policy/value 诊断图，不弹出交互窗口。
    pick_ages = unique(max(1, min(par.JJ, [1, 20, par.JR, 60, par.JJ])));
    ip = 1;
    is = ceil(par.NS/2);

    fig = figure('Visible', 'off');
    hold on
    for k = 1:numel(pick_ages)
        j = pick_ages(k);
        plot(col.a_nodes, sol(j).aplus(:, ip, is), 'DisplayName', sprintf('Age %d', ages(j)));
    end
    xlabel('Assets');
    ylabel('Next assets');
    legend('Location', 'best');
    title('Savings policy');
    saveas(fig, 'labor_collocation_policy_aplus.png');
    % close(fig);

    fig = figure('Visible', 'off');
    hold on
    for k = 1:numel(pick_ages)
        j = pick_ages(k);
        plot(col.a_nodes, sol(j).c(:, ip, is), 'DisplayName', sprintf('Age %d', ages(j)));
    end
    xlabel('Assets');
    ylabel('Consumption');
    legend('Location', 'best');
    title('Consumption policy');
    saveas(fig, 'labor_collocation_policy_consumption.png');
    % close(fig);

    fig = figure('Visible', 'off');
    hold on
    for k = 1:numel(pick_ages)
        j = pick_ages(k);
        plot(col.a_nodes, sol(j).l(:, ip, is), 'DisplayName', sprintf('Age %d', ages(j)));
    end
    xlabel('Assets');
    ylabel('Hours');
    legend('Location', 'best');
    title('Labor policy');
    saveas(fig, 'labor_collocation_policy_labor.png');
    % close(fig);

    fig = figure('Visible', 'off');
    hold on
    for k = 1:numel(pick_ages)
        j = pick_ages(k);
        plot(col.a_nodes, sol(j).V(:, ip, is), 'DisplayName', sprintf('Age %d', ages(j)));
    end
    xlabel('Assets');
    ylabel('Value');
    legend('Location', 'best');
    title('Value function');
    saveas(fig, 'labor_collocation_value.png');
    % close(fig);

    fig = figure('Visible', 'off');
    semilogy(ages, max(residual, eps), 'LineWidth', 1.25);
    xlabel('Age');
    ylabel('Infinity norm');
    title('Newton residual by age');
    saveas(fig, 'labor_collocation_residual.png');
    % close(fig);
end

function make_distribution_plots(agg, dist, ages)
    % 保存分布推进后 cohort profile 图。
    fig = figure('Visible', 'off');
    plot(ages, agg.c_coh, 'LineWidth', 1.25);
    hold on
    plot(ages, agg.y_coh + agg.pen_coh, 'LineWidth', 1.25);
    xlabel('Age');
    ylabel('Mean');
    legend({'Consumption', 'Earnings plus pension'}, 'Location', 'best');
    title('Cohort means');
    saveas(fig, 'labor_collocation_cohort_means.png');
    % close(fig);

    fig = figure('Visible', 'off');
    plot(ages, agg.l_coh, 'LineWidth', 1.25);
    xlabel('Age');
    ylabel('Hours');
    title('Cohort labor');
    saveas(fig, 'labor_collocation_cohort_labor.png');
    % close(fig);

    fig = figure('Visible', 'off');
    plot(ages, agg.a_coh, 'LineWidth', 1.25);
    xlabel('Age');
    ylabel('Assets');
    title('Cohort assets');
    saveas(fig, 'labor_collocation_cohort_assets.png');
    % close(fig);

    fig = figure('Visible', 'off');
    plot(ages, agg.cv_c, 'LineWidth', 1.25);
    hold on
    plot(ages, agg.cv_y, 'LineWidth', 1.25);
    xlabel('Age');
    ylabel('Coefficient of variation');
    legend({'Consumption', 'Earnings'}, 'Location', 'best');
    title('Cohort dispersion');
    saveas(fig, 'labor_collocation_cohort_cv.png');
    % close(fig);

    fig = figure('Visible', 'off');
    plot(ages, agg.cv_y, 'LineWidth', 1.25);
    hold on
    plot(ages, agg.cv_l, 'LineWidth', 1.25);
    plot(ages, agg.cv_h, 'LineWidth', 1.25);
    plot(ages, agg.corr_hl, 'LineWidth', 1.25);
    xlabel('Age');
    ylabel('Variance decomposition');
    legend({'Earnings', 'Hours', 'Productivity', 'Correlation'}, 'Location', 'best');
    title('Variance decomposition earnings');
    saveas(fig, 'labor_collocation_variance_decomposition.png');
    % close(fig);

    fig = figure('Visible', 'off');
    plot(ages, agg.frac_bor, 'LineWidth', 1.25);
    xlabel('Age');
    ylabel('Fraction');
    title('Borrowing constrained households');
    saveas(fig, 'labor_collocation_borrowing_constrained.png');
    % close(fig);

    fig = figure('Visible', 'off');
    semilogy(ages, max(dist.mass_error, eps), 'LineWidth', 1.25);
    xlabel('Age');
    ylabel('Absolute error');
    title('Distribution mass error');
    saveas(fig, 'labor_collocation_distribution_mass_error.png');
    % close(fig);
end
