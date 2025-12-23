module Data_pre
export Detrend, Common_Detrend, Difference_method, HPfilter
using  LinearAlgebra, Statistics, SparseArrays


function Detrend(data)
    n = length(data)
    t = collect(1:n)
    X = hcat(ones(n), t.-1)
    β = X \ data
    trend = X * β
    detrended_data = data - trend
    return detrended_data,trend
end



function Common_Detrend(m_Data)
    T,m = size(m_Data)
    tt = kron(ones(m), 0:T-1) 
    X = kron(I(m), ones(T))
    Y = vec(m_Data)
    Z = hcat(X,tt)

    θ = Z \ Y
    Y_hat = Z*θ
    ϵ_hat =  Y - Y_hat
    trend = reshape(Y_hat,T,m)
    detrend_data = reshape(ϵ_hat,T,m)
    return detrend_data, trend
end

function Difference_method(v_data)
    y = v_data[1:end-1]
    yp = v_data[2:end]
    Δ_y = yp .- y
    γ = mean(Δ_y)
    ϵ_hat = Δ_y .- γ
    return ϵ_hat, γ
end

function HPfilter(y; weight = 1600)
    # Hodrick-Prescott filter implementation
    # weight: smoothing parameter, default is 1600 for quarterly data, 129600 for monthly data,6.25 for annual data
     λ = float(weight)
    T = size(y, 1)
    @assert T ≥ 4 "HP filter needs at least 4 observations."

    d2 = fill(λ, T-2)                         # ±2 对角
    d1 = fill(-4λ, T-1); d1[1] = -2λ; d1[end] = -2λ   # ±1 对角
    d0 = fill(1 + 6λ, T)                      # 主对角
    d0[1] = 1 + λ
    d0[2] = 1 + 5λ
    d0[end-1] = 1 + 5λ
    d0[end] = 1 + λ

    B = spdiagm(-2 => d2, -1 => d1, 0 => d0, 1 => d1, 2 => d2)  # 对称五对角
    trend = B \ y
    cycle = y .- trend
    return cycle, trend
end

end
