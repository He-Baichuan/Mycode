using LinearAlgebra
using Interpolations
"""
distributionSS(P::AbstractMatrix{<:Real}; method=:eigen, maxit = 1e5)
基于特征向量法和模拟方法计算状态转移矩阵 P 的稳态分布 π，
要求 P 是行和为 1 的方阵。
"""
function distributionSS(P::AbstractMatrix{<:Real}; method=:eigen,maxit = 1e5)
    P = transpose(P)  #先转置为方便求解的矩阵形式
    n, m = size(P)
    @assert n == m "P 必须是方阵"
    # 检验转置后的矩阵列和为 1，也就是说原始状态转移矩阵的行和为 1
    tol = 1e-12
    @assert all(abs.(sum(P, dims=1)[:] .- 1) .< tol) "P 的每行和必须为 1"
    if method == :eigen
        # 对 P 求右特征分解
        F = eigen(P)
        # 定位特征值最接近 1 的索引
        idx = argmin(abs.(F.values .- 1))
        # 提取对应的特征向量，并取实部
        v = real(F.vectors[:, idx])
        # 归一化特征向量
        π = v / sum(v)
        return π
    elseif method == :simulate
        v = zeros(size(P, 1))
        v[1] = 1.0;
        for it in 1:maxit
            v = P * v
           if it % 100 ==0 && all(abs.(v-P*v).<1e-12)
                break
            end
        end
        return v

    else
        error("不支持的方法：$method")
    end
end

function interp(x,y,x1;right_extrap=false, left_extrap=false)
    interpval = linear_interpolation(x,y,extrapolation_bc = Line());
    y1 = interpval(x1);
    
    if ~right_extrap 
        above = x1 .> x[end];
        y1[above] .= y[end];
    end
    if ~left_extrap
        below = x1 .< x[1];
        y1[below] .= y[1];
    end

    return y1
end


function SolveEGM(c0,Xt,par,grid)
    Va = uPrime(par,c0);
    #loop until convergence
    # println("solving for policy rules")
    tol = 1e-10
    
    for it =1:10000
        Va1 = EGMStepBack(Va,Xt,par,grid)[1]

        if it % 50 == 0
            test = abs.(Va1-Va)/(abs.(Va) .+tol)
            # println("it = $(it), test = $(maximum(test))")
            if all(test  .< tol)
                break
            end
        end

        Va = Va1
    end
    return EGMStepBack(Va,Xt,par,grid)
end

