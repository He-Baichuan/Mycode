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