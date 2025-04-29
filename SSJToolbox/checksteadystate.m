function checksteadystate(m, f, tol)
        %   CHECKSTEADYSTATE 检查模型稳态在全局堆叠系统下是否收敛
        %   checkSteadyState(m, f) 使用默认容忍度 tol = 1e-6
        %   checkSteadyState(m, f, tol) 指定容忍度
        %
        %   输入：
        %     m   - ModelUtils 环境对象，包含 m.ss.initial, m.ss.terminal, m.T, m.nx
        %     f   - 残差函数句柄，调用时用法 y = f(m, X, E)
        %     tol - 收敛容忍度（可选，默认 1e-6）

            if nargin < 3 || isempty(tol)
               tol = 1e-6;
            end

        % 1. 检查初末稳态是否一致
            diffST = abs(m.ss.initial - m.ss.terminal);
            if any(diffST > tol)
                warning('初始稳态与末期稳态在容忍度 %g 下不一致，跳过稳态检查。', tol);
                return;
            end

            % 2. 构造全局稳态路径
            [Xss, Ess] = longsteadystate(m);

             % 3. 计算全局残差
            res = f(m, Xss, Ess);

             % 4. 检查残差是否小于容忍度
            idx = find(abs(res) > tol);
            if ~isempty(idx)
             % 将索引转为 [period, eqn] 格式
                periods = floor((idx-1) / m.nx) + 1;
                eqns    = mod(idx-1, m.nx) + 1;
                 fprintf('以下期次/方程的残差超过 tol=%g：', tol);
                for k = 1:length(idx)
                    fprintf('  period %d, equation %d: residual = %.3e', periods(k), eqns(k), res(idx(k)));
                end
                error('稳态检验失败：存在残差超出容忍度。');
            end

            fprintf('稳态检验通过：最大残差 = %.3e <= tol = %g', max(abs(res)), tol);
        end