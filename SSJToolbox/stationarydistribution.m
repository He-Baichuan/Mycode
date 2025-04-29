function D = stationarydistribution(T,method,maxit)
% stationarydistribution  求解列随机矩阵 T 的平稳分布
%   要求 sum(T,1)==1
    tol = 1e-12;
    if any(abs(sum(T,1) - 1) > tol)
        error('stationarydistribution: each column of T must sum to 1');
    end
    if method == 0
    %——用特征值法——
    [V, Dmat] = eig(T);
    % 找到最接近1的特征值对应的特征向量
    [~, idx] = min(abs(diag(Dmat) - 1));
    v = real(V(:,idx));
    D = v / sum(v);
    elseif method == 1
    %----用模拟法----%
         n = size(T,1);
            D = zeros(n,1);
            D(1) = 1;  % 初始分布
            for it = 1:maxit
                Dnew = T * D;
                if mod(it,100)==0 && max(abs(Dnew - D)) < 1e-8
                    D = Dnew;
                    break
                end
                D = Dnew;
            end
            D = D / sum(D);

    else
            error('未知方法: %s. 支持 ''eigen'' 或 ''simulate''.', method);
    end
end