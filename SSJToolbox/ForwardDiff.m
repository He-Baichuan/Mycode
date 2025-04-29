function [FX, FE] = ForwardDiff(f, X, E, T, varargin)
    % if nargin < 10
        h = 1e-6;
    % end
    n = length(X);
    m = length(E);
    f0 = f(X, E, T,varargin);
    FX = zeros(n, n);
    FE = zeros(n, m);

    for i = 1:n
        X_perturbed = X;
        X_perturbed(i) = X_perturbed(i) + h;
        FX(:, i) = (f(X_perturbed, E, T, varargin) - f0) / h;
    end
    for j = 1:m
        E_perturbed = E;
        E_perturbed(j) = E_perturbed(j) + h;
        FE(:, j) = (f(X, E_perturbed, T, varargin) - f0) / h;
    end
    if rcond(FX) < 1e-12
    error('矩阵 F_X 奇异，无法求解线性系统。请检查模型设定。');
    end
end
