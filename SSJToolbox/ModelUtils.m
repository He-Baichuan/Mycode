% ModelUtils.m: MATLAB translation of the Julia `ModelUtils` module
% Integrates full functionality: model environment setup, steady-state handling,
% indexing, reshaping, lag/lead operations, Jacobians, IRFs, nonlinear and optimal policy.

classdef ModelUtils
    properties
        par            % user-defined parameters
        vars           % containers.Map: keys {'0','-1','1','exog',-2,2,...} -> cell array of var names
        T              % time horizon
        ss             % struct with fields: initial (nx×1), terminal (nx×1), exog (nexog×1)
        nx             % number of endogenous variables
        nexog          % number of exogenous shocks
        ss_sym         % symbolic steady-state (for non-linear policy)
        idx_now         % >>> 新增：当期内生变量的整数索引（cell->double）
        idx_exog        % >>> 新增：外生变量索引
    end

    methods (Static)
        function m = ModelEnv(par, vars, steadystate, T)
            % Constructor: initialize ModelUtils environment
            % varFields: struct with fields '0','-1','1',...,'exog' mapping to cell arrays of names
            m = ModelUtils;
            m.par = par;
            % m.vars = containers.Map();
            m.vars.endo = vars.endo;
            m.vars.exog = vars.exog;
            % for k = fieldnames(vars)'
            %     m.vars(k{1}) = vars.(k{1});
            % end
            m.T = T;
            m.nx = numel(steadystate.initial);
            m.nexog = numel(steadystate.exog);
            % ensure terminal steady state
            if ~isfield(steadystate,'terminal')
                steadystate.terminal = steadystate.initial;
            end
            m.ss = steadystate;
            % m.ss_sym = struct();
            m.idx_now  = containers.Map(vars.endo,  1:m.nx);
            m.idx_exog = containers.Map(vars.exog,1:m.nexog);
        end

        % function val = steadystatevalue(s, m, endogexog, which)
        %     % Return steady-state value of variable s
        %     if nargin<4, which = 'terminal'; end
        %     if strcmp(endogexog,'endog')
        %         vec = m.ss.(which);
        %         names = m.vars('0');
        %     elseif strcmp(endogexog,'exog')
        %         vec = m.ss.exog;
        %         names = m.vars('exog');
        %     else
        %         error('endogexog must be ''endog'' or ''exog''.');
        %     end
        %     idx = find(strcmp(names, s)); assert(~isempty(idx),'Variable not found');
        %     val = vec(idx);
        % end

        % function idx = varindex(s, m, endogexog)
        %     % Return indices in long X/E for variable s
        %     if nargin<3, endogexog='endog'; end
        %     if strcmp(endogexog,'endog')
        %         names = m.vars('0'); n = m.nx;
        %     else
        %         names = m.vars('exog'); n = m.nexog;
        %     end
        %     j = find(strcmp(names, s)); assert(~isempty(j),'Variable not found');
        %     idx = (j-1)*m.T + (1:m.T);
        % end

        % function checksteadystate(m, f, tol)
        %     % Check f(m,X_ss,E_ss)=0 within tolerance
        %     if nargin<3
        %         tol=1e-9;
        %     end
        %     if any(abs(m.ss.initial-m.ss.terminal)>tol)
        %         warning('Initial and terminal steady states differ; skipping check.');
        %         return;
        %     end
        %     [Xss, Ess] = ModelUtils.longsteadystate(m);
        %     res = f(m,Xss,Ess);
        %     if any(abs(res)>tol)
        %         error('Steady-state residuals exceed tol.');
        %     end
        % end

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
            [Xss, Ess] = ModelUtils.longsteadystate(m);

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


        function [Xss, Ess] = longsteadystate(m)
            % Long-format steady state: repeat init/terminal and exog across T
            Xss = repmat(m.ss.terminal(:), m.T, 1);
            Ess = repmat(m.ss.exog(:), m.T, 1);
        end

        function blocks = wide(x, n)
            % Split long vector x (n*T×1) into cell array of n vectors length T
            T = numel(x)/n;
            x = reshape(x,n,T);
            x = x';
            blocks = mat2cell(x,T, ones(1,n));
        end

        % function xlong = long(blocks)
        %     % Concatenate cell array of vectors into long vector
        %     xlong = vertcat(blocks{:});
        % end
        % 
        function S = id(name,m,num)
            if num == 0
                S = m.idx_now(name);
            % endo = m.vars('0');
            % S = find(strcmp(endo,name));
            end
            if num == -1
                S= m.idx_now([name '_l']);
                % lag = m.vars('-1');
                % S = find(strcmp(lag,name));
            end
            if num ==1
                S = m.idx_now([name '_p']);
                % lead = m.vars('-1');
                % S = find(strcmp(lead,name));
            end
        end
    
        % function id = id(name, m,leadlag)
        %     % leadlag =  0 当期，<0 滞后，>0 前导
        %     switch leadlag
        %         case 0,  id = m.idx_now(name);
        %         case -1, id = m.idx_now([name '_l']);     % 依据命名规则
        %         case 1,  id = m.idx_now([name '_p']);
        %         otherwise, error('暂仅支持 ±1 期');
        %     end
        % end
    
        % function S = contemp(X, m)
        %     % Return struct of contemporaneous values
        %     names = m.vars('0'); 
        %     blocks = ModelUtils.wide(X,m.nx);
        %     for i=1:m.nx
        %         S.(names{i}) = blocks{i}(1); 
        %     end
        % end
        % 
         function S = contemp(X,name,m)
            idx =  ModelUtils.id(name,m,0);
            S = X(idx:m.nx:end);
         end
        % function S = lag(X, m, j)
        %     % Return struct of lagged values (default j=1)
        %     if nargin<3, j=1; end
        %     names = m.vars(['-' num2str(j)]); % require varFields has key '-j'
        %     blocks = ModelUtils.wide(X,m.nx);
        %     % shift each by prepending steady-state initial
        %     for i=1:m.nx
        %         xi = blocks{i};
        %         S.(names{i}) = [repmat(m.ss.initial(i),j,1); xi(1:end-j)];
        %     end
        % end
        % 
        % function S_l = lag(var)
        %  S_l = lagmatrix(var,1);
        %  S_l(1) = var(1);
        % 
        % end
        % function S = lead(X, m, j)
        %     % Return struct of lead values (default j=1)
        %     if nargin<3, j=1; end
        %     names = m.vars(num2str(j)); % require varFields key 'j'
        %     blocks = ModelUtils.wide(X,m.nx);
        %     for i=1:m.nx
        %         xi = blocks{i};
        %         S.(names{i}) = [xi(j+1:end); repmat(m.ss.terminal(i),j,1)];
        %     end
        % end
        % function S_p = lead(var,T)
        %  S_p = lagmatrix(var,-1);
        %  S_p(T) = var(T);
        % 
        % end
        function v_l = lag(v)
            v_l      = circshift(v,1);
            v_l(1) = v(1);        % 第一行置稳态
        end
        function v_f = lead(v)
            v_f        = circshift(v,-1);
            v_f(end) =v(end);      % 末行置稳态
        end
    
        % function S = exogenous(E, m)
        %     % Return struct of exogenous series
        %     names = m.vars('exog'); blocks = ModelUtils.wide(E,m.nexog);
        %     for i=1:m.nexog, S.(names{i}) = blocks{i}(1); end
        % end

        function [FX, FE] = get_jacobian(f, m, X, E)

            f0 = f(m, X, E);
            n  = numel(X); 
            p = numel(E); 
            h = sqrt(eps);
            FX = zeros(numel(f0), n); 
            FE = zeros(numel(f0), p);
            for i = 1:n
                dX = zeros(n,1); 
                dX(i)=h;
                FX(:,i) = (f(m, X+dX, E) - f0)/h;
            end
            for j = 1:p
                dE = zeros(p,1); 
                dE(j)=h;
                FE(:,j) = (f(m, X, E+dE) - f0)/h;
            end
            if rcond(FX) < 1e-12
               error('矩阵 F_X 奇异，无法求解线性系统。请检查模型设定。');
            end
        end
        % function [FX,FE] = get_jacobian(f,m,X,E)
        %     f0 = f(m,X,E);            n = numel(X); p = numel(E);
        %     hX = ModelUtils.FD_STEP * max(1,abs(X));     % ❷ 比例步长
        %     hE = ModelUtils.FD_STEP * max(1,abs(E));
        % 
        %     % 一次性生成扰动矩阵，避免 n+p 次函数调用各自分配
        %     DX = spdiags(hX,0,n,n);                     % 稀疏对角
        %     DE = spdiags(hE,0,p,p);
        % 
        %     FX = (f(m, X+DX, E) - f0) ./ hX.';          % 隐式广播
        %     FE = (f(m, X, E+DE) - f0) ./ hE.';
        % 
        %     % 稳健奇异性检查
        %     if rcond(FX) < 1e-10
        %         error('F_X 奇异或病态 (rcond=%g)', rcond(FX));
        %     end
        % end


         function FX = get_fX(f, m, X, E)
            if nargin<4, [X,E] = ModelUtils.longsteadystate(m); end
            [FX, ~] = ModelUtils.get_jacobian(f, m, X, E);
        end

        function FE = get_fE(f, m, X, E)
            if nargin<4, [X,E] = ModelUtils.longsteadystate(m); end
            [~, FE] = ModelUtils.get_jacobian(f, m, X, E);
        end

        function IRF = linearIRFs(f, m, X, E)
            if nargin<3, [X,E] = ModelUtils.longsteadystate(m); end
            [FX, FE] = ModelUtils.get_jacobian(f, m, X, E);
            IRF =-( sparse(FX) \ FE );
        end

        function X = nonlineartransition(m, f, X0, E, maxit)
            % Solve non-linear transition via Newton iterations
            if nargin<5, maxit=100; end
            X = X0;
            for it=1:maxit
                Xold = X;
                res = f(m,X,E);
                FX = ModelUtils.get_fX(f,m,X,E);
                X = X - (FX \ res);
                update = max(abs(X-Xold));
                if update<1e-6, return; end
            end
            error('Nonlinear transition did not converge');
        end

        function plotIRFs(IRFMat, m, irf_horizon, shock, shock_horizon, series)
            % Plot IRFs from IRFMat (nx*T×nexog*T)
            if nargin<3
                irf_horizon=80; 
            end
            if nargin<4 
                shock=1;
            end
            if nargin<5
                shock_horizon=0; 
            end
            if nargin<6
                series=m.vars.endo; 
            end
            % determine shock index
            if isnumeric(shock)
                sidx=shock;
            else
                exn=m.vars.exog; 
                sidx=find(strcmp(exn,shock)); 
            end
            vec = IRFMat(:,(sidx-1)*m.T + shock_horizon+1);
            blocks = ModelUtils.wide(vec,m.nx);
            figure;
            for i=1:m.nx
                subplot(ceil(m.nx/3),3,i);
                plot(1:irf_horizon, blocks{i}(1:irf_horizon),'LineWidth', 2);
                xlabel('时间');
                title(series{i});
                sgtitle('脉冲响应函数');
            end
        end

        function [Xopt, lam] = optimaltransitionpath(m, Ugrad, Uhess, f, X0, E, maxit)
            % Compute optimal transition with policy via first-order conditions
            if nargin<7, maxit=100; end
            X = X0;
            FX = ModelUtils.get_fX(f,m,X,E);
            UX = Ugrad(m,X,E);
            lam = FX' \ UX;
            neq = size(FX,1);
            for it=1:maxit
                FX = ModelUtils.get_fX(f,m,X,E);
                UX = Ugrad(m,X,E);
                res1 = UX - FX'*lam;
                res2 = f(m,X,E);
                % approximate Hessian of Lagrangian via numeric
                H = ModelUtils.numericHessian(@(x) Ugrad(m,x,E) - FX'*lam, X);
                % build KKT Jacobian
                K = [H, -FX'; FX, zeros(neq)];
                upd = -K \ [res1; res2];
                X = X + upd(1:numel(X)); lam = lam + upd(numel(X)+1:end);
                if max(abs(upd)) < 1e-7, break; end
            end
            Xopt = X;
        end

        function [xopt, lam] = optimalLQpolicy(Q, R, A, B)
            % Exact LQ policy solver: min 0.5 x'Qx + R'x s.t. A + Bx = 0
            nX = size(Q,1);
            K = [Q, -B'; B, zeros(size(B,1))];
            sol = -K \ [R;A];
            xopt = sol(1:nX);
            lam  = sol(nX+1:end);
        end

        function [Xnew, lam] = optimalLQpolicy_point(m, Ugrad, Uhess, f, Ehat, approx)
            % LQ around approx point or steady state
            if nargin<6
                [X0, E0] = ModelUtils.longsteadystate(m);
            else
                X0 = approx{:}; E0 = approx{2};
            end
            Q = Uhess(m,X0,E0);
            R = Ugrad(m,X0,E0);
            B = ModelUtils.get_fX(f,m,X0,E0);
            A = ModelUtils.get_fE(f,m,X0,E0)*Ehat + f(m,X0,E0);
            [dx, lam] = ModelUtils.optimalLQpolicy(Q,R,A,B);
            Xnew = X0 + dx;
        end

        function H = numericHessian(fun, x)
            % Numerically approximate Hessian of vector function fun(x)'*lam is embedded above
            % Here fun returns gradient; Hessian approximated per element
            n = numel(x); h=1e-5;
            g0 = fun(x);
            H = zeros(n);
            for i=1:n
                ei = zeros(n,1); ei(i)=h;
                g1 = fun(x+ei);
                H(:,i) = (g1 - g0)/h;
            end
        end
    end
end
