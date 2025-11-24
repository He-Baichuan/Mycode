using Plots
using Parameters
using LinearAlgebra
using LinearInterpolations
using Roots
using Random
using Distributions
include("Numerical_method.jl")
#首先我们设定模型的参数
function SettingPar_TransitionDynamics(;
    nA = 300,
    nS = 7,
    α = 0.33,
    β = 0.96,
    γ = 2.0,
    δ = 0.05,
    ϕ = 0.0,
    ρ = 0.9,
    σ = 0.03,
    Zbar = 0.0, # 总冲击log
    ρ_z = 0.95,
    T = 150,
    # Tax params.
    ShockSize = 0.01,
    τ_l0 = 0.0,
    τ_l1 = 0.0,
    τ_k0 = 0.0,
    τ_k1 = 0.0,
    G0 = 0.0,
    G1 = 0.0,
    a_max = 250.0,
    ν = 0.025)

    """
    生成劳动格点
    """
    gridS, transS = MyNumerical_method.rouwenhorst(nS, ρ, σ)
    gridS = exp.(gridS)

    """
    计算不变分布
    """
    inv_S = ones(nS)/nS#初始化
    for i in 1:10000 
        inv_S_new = (inv_S' * transS)'
        if maximum(abs.(inv_S_new - inv_S))<1e-8
            break
        end
        inv_S .= inv_S_new
    end
    Lbar = sum(gridS.*inv_S)#劳动禀赋
    
    """
    生成资产格点
    """
    if ν == 0.0
        gridA = collect(range(-ϕ, a_max, nA))
    else
        step = 0:(nA - 1)
        gridA = -ϕ .+ (a_max + ϕ).*((1.0 + ν).^(step) .- 1.0)./((1.0 + ν)^(nA - 1) - 1.0)
    end

    """
    构造冲击序列
    """
    Shock_Z = zeros(T)
    Shock_Z[1] = Zbar; Shock_Z[T] = Zbar; Shock_Z[2] = ShockSize
    for t in 3:T-1
        Shock_Z[t] = ρ_z * Shock_Z[t-1]
    end
    Shock_Z = exp.(Shock_Z)
    Zbar = exp(Zbar)

    """
    构造政策冲击序列（永久性变化）
    """
    τl_seq = zeros(T); τl_seq .= τ_l1; τl_seq[1] = τ_l0
    τk_seq = zeros(T); τk_seq .= τ_k1; τk_seq[1] = τ_k0
    G_seq = zeros(T); G_seq .= G1; G_seq[1] = G0
return(α = α, β = β, γ = γ, δ = δ, ϕ = ϕ, Zbar = Zbar, ρ_z = ρ_z,
        gridA = gridA, gridS = gridS, transS = transS, Lbar = Lbar,
        T = T, nA = nA, nS = nS, Shock_Z = Shock_Z, τl_seq = τl_seq,
        τk_seq = τk_seq, Gseq = G_seq)
end

function SolveHHprobEGM(param, r, w, τl, τk, transfer)
    @unpack β, δ, γ, nA, nS, gridA, gridS, transS = param 
    gAA = repeat(gridA, 1, nS)
    gSS = repeat(gridS', nA, 1)
    gYY = (1.0 + r*(1.0 -τk)).*gAA .+ w*(1.0 -τl).*gSS .+ transfer# cash on hand，可支配收入
    tol = 1e-8
    maxit = 5000
    cp = (gYY - gAA)/2
    ap = copy(cp)
    A_endo = copy(cp)
    function EGMIteration(cp)
        RHS = β*(1.0 +r*(1-τk))* cp.^(-γ)*transS'
        c = RHS.^(-1.0/γ)
        A_endo = (c + gAA - w*(1-τl).*gSS .- transfer)/(1.0 + r*(1.0 - τk))
        ap_new = zeros(nA, nS)
        for j in 1:nS
            ap_new[:,j] = MyNumerical_method.linearGridIntp!(A_endo[:,j], gridA, gridA, ap_new[:,j])
        
            for i in 1:nA
                if ap_new[i,j]<gridA[1]
                    ap_new[i,j] = gridA[1]
                else
                    break
                end
            end
        end
        cp_new = gYY - ap_new
        return cp_new, ap_new, A_endo
    end
    for iter in 1:maxit
        cp_new, ap, A_endo = EGMIteration(cp)
        crit = maximum(abs.(cp_new - cp))
        if crit<tol
            println("HH Problem is solved by EGM ")
            break
        end
        if iter == maxit
            error("The max iteration !")
        end
        cp .= cp_new
    end
    return (cp = cp, ap = ap, A_endo = A_endo)
end

function SolveInvariantDist(param, decision)
    @unpack transS, nA, nS, gridA = param 
    @unpack ap = decision
    maxit = 500000
    tol = 1e-8
    id_L = zeros(Int,nA,nS)
    weight_L = zeros(nA,nS)
    for j in 1:nS
        id_L[:,j], weight_L[:,j] = MyNumerical_method.coordGridIntp!(gridA, ap[:,j],id_L[:,j], weight_L[:,j], robust = true)
    end

    dsn = fill(1.0/(nA*nS),(nA, nS))
    for iter in 1:maxit
        dsn_new = zeros(nA, nS)
        for i in 1:nA, j in 1:nS #构造内层概率
            if dsn[i,j]>0.0
                dsn_new[id_L[i,j],j] += dsn[i,j]*weight_L[i,j]
                dsn_new[id_L[i,j]+1,j] += dsn[i,j]*(1.0 - weight_L[i,j])
            end
        end
        dsn_new = dsn_new*transS
        crit = maximum(abs.(dsn_new - dsn))
        dsn .= dsn_new
        if crit <tol
            println("Invarant Dist. is solved")
            break
        end
        if iter == maxit
            error("The Invarant Dist. is not converged")
        end
    end
    return dsn
end
function ExcessDemand(param, VaratTime, Kguess)
    @unpack Lbar, α, δ, nA, nS, gridA = param
    Zbar = VaratTime[1]
    τl= VaratTime[2]
    τk = VaratTime[3]
    Gexp = VaratTime[4]

    r = Zbar*α *Kguess^(α-1)*Lbar^(1-α)-δ
    w = Zbar*(1-α)*Kguess^α*Lbar^(-α)
    Y = Zbar*Kguess^α*Lbar^(1-α)
    G = Y*Gexp
    transfer = τl*w*Lbar + τk*r*Kguess - G
    decision = SolveHHprobEGM(param, r, w, τl, τk, transfer)
    dsn = SolveInvariantDist(param, decision)
    gAA = repeat(gridA,1,nS)
    Ea = sum(dsn.*gAA)
    #excess = 2*(Ea - Kguess)/(Ea + Kguess)
    excess = Ea - Kguess
    return (excess, dsn, decision, w,r , transfer, Ea) 
end

function SolveModel(param, VaratTime)
    @unpack α, δ, β, Lbar = param 
    Zbar = VaratTime[1]
    Kcomplete = (α/(1/β -1 +δ))^(1/(1-α))*Lbar
    r_low = 0.001
    Kmax = Zbar*(α/(r_low+δ))^(1/(1-α))*Lbar
    	function objFct(Kguess);
		println("\nCapital Guess: $Kguess")
		(excDem, ) = ExcessDemand(param, VaratTime, Kguess);
		println("\nExcess Demand: $excDem")
		return excDem
	    end
    eqm_K = fzero(k -> objFct(k), Kcomplete, Kmax, Roots.Brent(),atol = 0.001)
    #eqm_K = MyNumerical_method.Rootfinding_Bisection(objFct, 0.1, 2*Kmax)
    (excess, dsn, decision, w, r, transfer, Ea)  = ExcessDemand(param, VaratTime, eqm_K)
    return (decision, dsn,   w, r, transfer,eqm_K, Ea)
end

"""
价格方程：给定 (K,Z) 计算 r,w
Y = Z K^α L^{1-α}, r = α Y/K - δ, w = (1-α) Y/L
"""
function prices_Aiy(K::Real, Z::Real, param)
    @unpack α, δ, Lbar = param
    Y = Z * K^α * Lbar^(1-α)
    r = α * Y / K - δ
    w = (1-α) * Y / Lbar
    return r, w
end

"""
财政方程：给定 (K,Z,r,w)，用参数中的 τ_l, τ_k, G/Y 计算总转移 transfer
transfer = τ_l w L + τ_k r K - G
"""
function fiscal_Aiy(K::Real, Z::Real, r::Real, w::Real, param)
    @unpack Lbar, α, τl_seq, τk_seq, Gseq = param
    τl   = τl_seq[1]
    τk   = τk_seq[1]
    Gexp = Gseq[1]          # G/Y 的稳态值
    Y    = Z * K^α * Lbar^(1-α)
    G    = Gexp * Y
    transfer = τl*w*Lbar + τk*r*K - G
    return τl, τk, transfer
end

"""
一步分布转移算子: 从 (D_L, policy ap) 推出下一期分布 D_t
基本上就是把你 SolveDistSeq 里的单期循环抽出来
"""
function OneStepDist(param, ap::AbstractMatrix, D_L::AbstractMatrix)
    @unpack nA, nS, gridA, transS = param

    # ap 里可能有 Dual，D_L 是当前分布（Float 或 Dual）
    # 统一用一个“提升后”的数值类型，确保能存下 Dual
    Tprom = promote_type(eltype(ap), eltype(D_L))

    id_L     = zeros(Int,   nA, nS)
    weight_L = zeros(Tprom, nA, nS)

    # 这里用的是你在 Numerical_method.jl 里定义的 coordGridIntp_generic!
    for j in 1:nS
        MyNumerical_method.coordGridIntp_generic!(
            gridA,
            view(ap, :, j),        # X，可以是 Dual
            view(id_L, :, j),
            view(weight_L, :, j)   # 权重是 Tprom（能容纳 Dual）
        )
    end

    dsn_new = zeros(Tprom, nA, nS)
    for i in 1:nA, j in 1:nS
        mass = D_L[i,j]           # 可能是 Float 或 Dual
        if mass != zero(mass)
            dsn_new[id_L[i,j],   j] += mass * weight_L[i,j]
            dsn_new[id_L[i,j]+1, j] += mass * (one(Tprom) - weight_L[i,j])
        end
    end

    # 乘以 transS（Float64），返回类型仍然是 Tprom
    D_t = dsn_new * transS
    return D_t
end

"""
Euler 方程残差：
给定当前/下一期政策 g_t,g_{t+1} (都是保存为 ap = a'(a,s) 的矩阵) 和 (K,Z,K',Z'),
构造在每个 (a_i,s_j) 上
u_c(c_t) - β(1+r_{t+1}(1-τ_k)) E[u_c(c_{t+1})] = 0 的残差。
"""
function EulerResidual_Aiy(g_vec::AbstractVector,
                           gP_vec::AbstractVector,
                           K::Real, Z::Real,
                           KP::Real, ZP::Real,
                           param)

    @unpack nA, nS, gridA, gridS, transS, β, γ = param

    ap  = reshape(g_vec,  nA, nS)
    apP = reshape(gP_vec, nA, nS)

    gAA = repeat(gridA, 1, nS)
    gSS = repeat(gridS', nA, 1)

    # 当前价格与转移
    r, w = prices_Aiy(K, Z, param)
    τl, τk, trans = fiscal_Aiy(K, Z, r, w, param)

    # 下一期价格与转移
    rP, wP = prices_Aiy(KP, ZP, param)
    τlP, τkP, transP = fiscal_Aiy(KP, ZP, rP, wP, param)

    # 当前和下一期消费（政策函数定义在 “当前资产网格” 上）
    c  = (1 .+ r*(1-τk)).*gAA .+ w*(1-τl).*gSS .+ trans  .- ap
    cP = (1 .+ rP*(1-τkP)).*gAA .+ wP*(1-τlP).*gSS .+ transP .- apP

    # 防止数值误差导致负消费
    c  = max.(c,  1e-12)
    cP = max.(cP, 1e-12)

    u_c  = c.^(-γ)
    u_cP = cP.^(-γ)

    RHS = β * (1 .+ rP*(1-τkP)) .* (u_cP * transS')  # 对 s' 取期望

    res = u_c .- RHS
    return reshape(res, nA*nS)
end

"""
分布残差：
D_t - T(g_t) D_{t-1} = 0
X 里 D 用向量保存，这里 reshape 成矩阵
"""
function DistResidual_Aiy(g_vec::AbstractVector,
                          D_L_vec::AbstractVector,
                          D_vec::AbstractVector,
                          param)

    @unpack nA, nS = param
    ap   = reshape(g_vec,   nA, nS)
    D_L  = reshape(D_L_vec, nA, nS)
    D    = reshape(D_vec,   nA, nS)

    D_next = OneStepDist(param, ap, D_L)
    res = D .- D_next
    return reshape(res, nA*nS)
end

"""
总体残差：
1) 资本一致性 K_t = ∑ a_i D_t(a_i,s_j)
2) TFP AR(1): log Z_t - ρ_z log Z_{t-1} - ε_t = 0
(这里把 shock 的标准差规范为 1)
"""
function AggResidual_Aiy(D_vec::AbstractVector,
                         K::Real, Z::Real,
                         D_L_vec::AbstractVector,
                         Z_L::Real,
                         epsZ::Real,
                         param)

    @unpack nA, nS, gridA, ρ_z = param

    D   = reshape(D_vec,   nA, nS)
    D_L = reshape(D_L_vec, nA, nS)

    gAA = repeat(gridA, 1, nS)
    K_from_D = sum(D .* gAA)

    K_res = K - K_from_D
    Z_res = log(Z) - ρ_z*log(Z_L) - epsZ

    return [K_res; Z_res]
end

"""
利用你现有的 SolveModel 得到 Aiyagari 模型的稳态
返回稳态状态向量 xss = [vec(ap_ss); vec(D_ss); K_ss; Z_ss]
"""
function SteadyState_Aiy(param)
    # 使用设置中的稳态值 Zbar, τl_seq[1], τk_seq[1], Gseq[1]
    Var = [param.Zbar, param.τl_seq[1], param.τk_seq[1], param.Gseq[1]]
    decision, dsn, w_ss, r_ss, transfer_ss, K_ss, Ea_ss = SolveModel(param, Var)
    Z_ss = param.Zbar

    @unpack nA, nS = param
    g_ss = vec(decision.ap)
    D_ss = vec(dsn)

    xss = vcat(g_ss, D_ss, K_ss, Z_ss)
    return xss, decision, dsn, K_ss, Z_ss
end

"""
大系统 F_Aiy(X_L, X, X_P, ε) = 0
X = [vec(ap); vec(D); K; Z]
ε 目前是一维 TFP 冲击
"""
function F_Aiy(X_L::AbstractVector,
               X::AbstractVector,
               X_P::AbstractVector,
               epsilon::AbstractVector,
               param)

    @unpack nA, nS = param
    m  = nA*nS          # policy 维度
    md = m              # 分布维度

    # 拆 X_{t-1}
    g_L  = X_L[1:m]
    D_L  = X_L[m+1:m+md]
    K_L  = X_L[m+md+1]
    Z_L  = X_L[m+md+2]

    # 拆 X_t
    g    = X[1:m]
    D    = X[m+1:m+md]
    K    = X[m+md+1]
    Z    = X[m+md+2]

    # 拆 X_{t+1}
    g_P  = X_P[1:m]
    D_P  = X_P[m+1:m+md]
    K_P  = X_P[m+md+1]
    Z_P  = X_P[m+md+2]

    epsZ = epsilon[1]

    euler_root = EulerResidual_Aiy(g, g_P, K, Z, K_P, Z_P, param)
    dist_root  = DistResidual_Aiy(g, D_L, D, param)
    agg_root   = AggResidual_Aiy(D, K, Z, D_L, Z_L, epsZ, param)

    return vcat(euler_root, dist_root, agg_root)
end

"""
Rendahl 风格的 time-iteration 解线性系统:
A x_{t+1} + B x_t + C x_{t-1} + E ε_{t+1} = 0
求出 x_{t+1} = P x_t + Q ε_{t+1}
"""
function SolveSystem(A,B,C,E; maxit = 1000, tol = 1e-10)
    P = zeros(size(A))
    S = zeros(size(A))
    for it in 1:maxit
        P_new = -(A*P + B) \ C
        S_new = -(C*S + B) \ A
        res   = maximum(abs.(C + B*P_new + A*P_new*P_new))
        @info "Time iteration step $it, max residual = $res"
        P, S = P_new, S_new
        res < tol && break
    end

    Q = -(A*P + B) \ E

    # 简单稳定性检查
    eigP = eigen(P).values
    eigS = eigen(S).values
    if maximum(abs.(eigP)) > 1.0
        error("Non-existence: max |eig(P)| = $(maximum(abs.(eigP)))")
    end
    if maximum(abs.(eigS)) > 1.0
        error("No stable equilibrium: max |eig(S)| = $(maximum(abs.(eigS)))")
    end

    return P,Q
end

"""
主函数：运行 Reiter 方法
返回 (P,Q,xss,param,K_ss,Z_ss)
"""
function run_reiter_aiy()
    # 只需要稳态参数，ShockSize 对稳态无影响
    param = SettingPar_TransitionDynamics()

    xss, decision_ss, dsn_ss, K_ss, Z_ss = SteadyState_Aiy(param)

    # 在线性化点上计算 A,B,C,E
    A = ForwardDiff.jacobian(t -> F_Aiy(xss, xss, t, [0.0], param), xss)
    B = ForwardDiff.jacobian(t -> F_Aiy(xss, t, xss, [0.0], param), xss)
    C = ForwardDiff.jacobian(t -> F_Aiy(t, xss, xss, [0.0], param), xss)
    E = ForwardDiff.jacobian(t -> F_Aiy(xss, xss, xss, t, param), [0.0])
    #P, Q = t_iteration(C,B,A)
    #Q = Q*E
    #P, Q = SolveSystem(A,B,C,E)
    G,H,E,EE = TurnABCEtoSims(A,B,C,E)
    @time eu,G1,Impact = SolveQZ(G,H,E,EE)
    return ( P = G1, Q = Impact, xss = xss, param = param,
            K_ss = K_ss, Z_ss = Z_ss)
end


"""
画一个简单的 TFP 冲击 IRF（只看 K 和 Z）
"""
#=
function irf_aiy(result; T = 40)
    P    = result.P
    Q    = result.Q
    xss  = result.xss
    param = result.param
    K_ss = result.K_ss
    Z_ss = result.Z_ss

    @unpack nA, nS = param
    m  = nA*nS
    md = m
    idxK = 2*m + 1
    idxZ = 2*m + 2

    nx  = size(P,1)              # 注意：包含 TurnABCEtoSims 加的辅助状态
    IRF = zeros(nx, T)
    IRF[:,1] = Q[:,1]            # 一期 1 个单位冲击

    for t in 2:T
        IRF[:,t] = P * IRF[:,t-1]
    end

    K_irf = IRF[idxK,:] ./ K_ss .* 100.0
    Z_irf = IRF[idxZ,:] ./ Z_ss .* 100.0

    tgrid = 1:T
    p1 = plot(tgrid, Z_irf, lw=2, xlabel="t", ylabel="% dev from SS", label="Z")
    p2 = plot(tgrid, K_irf, lw=2, xlabel="t", ylabel="% dev from SS", label="K")
    plot(p1,p2,layout=(2,1),size=(800,600))
end
=#
function irf_aiy_times(result; T = 40)
    P    = result.P
    Q    = result.Q
    xss  = result.xss
    param = result.param
    K_ss = result.K_ss
    Z_ss = result.Z_ss

    @unpack nA, nS = param
    m  = nA*nS
    md = m
    idxK = 2*m + 1
    idxZ = 2*m + 2

    nx  = length(xss)
    IRF = zeros(nx, T)
    IRF[:,1] = Q[:,1]          # ε_1 = 1, 其他期 0

    for t in 2:T
        IRF[:,t] = P*IRF[:,t-1]
    end

    K_irf = IRF[idxK,:] ./ K_ss .* 100.0
    Z_irf = IRF[idxZ,:] ./ Z_ss .* 100.0

    tgrid = 1:T
    p1 = plot(tgrid, Z_irf, lw=2, xlabel="t", ylabel="% dev from SS", label="Z")
    p2 = plot(tgrid, K_irf, lw=2, xlabel="t", ylabel="% dev from SS", label="K")
    plot(p1,p2,layout=(2,1),size=(800,600))
end
function coordGridIntp_generic!(Grid::AbstractVector,
                                    X::AbstractVector,
                                    id_L::AbstractVector{Int},
                                    weight_L::AbstractVector)
        M = length(Grid)
        @assert M ≥ 2
        @assert length(X) == length(id_L) == length(weight_L)

        T = eltype(weight_L)   # e.g. Float64 或 Dual

        for i in eachindex(X)
            x = X[i]

            if x ≤ Grid[1]
                id_L[i] = 1
                weight_L[i] = one(T)          # 全部质量给左端点
            elseif x ≥ Grid[end]
                id_L[i] = M - 1
                weight_L[i] = zero(T)         # 全部质量给右端点
            else
                # 找到 k 使 Grid[k] ≤ x ≤ Grid[k+1]
                k = searchsortedlast(Grid, x)
                id_L[i] = k
                gL = Grid[k]
                gR = Grid[k+1]
                # 线性插值：左端点权重
                weight_L[i] = (gR - x) / (gR - gL)
            end
        end
        return id_L, weight_L
    end

    function t_iteration(A::AbstractMatrix,
                     B::AbstractMatrix,
                     C::AbstractMatrix;
                     μ::Real = 1e-6,
                     tol::Real = 1e-10,
                     maxit::Int = 10000)

    n = size(A, 1)
    T = promote_type(eltype(A), eltype(B), eltype(C))
    Iμ = μ * Matrix{T}(I, n, n)

    Ch = C
    Bh = B + C * (2Iμ)
    Ah = C * (Iμ * Iμ) + B * Iμ + A

    # 可选：奇异性检查
    if cond(Ah) > 1e16
        @warn "Matrix Ah is close to singular, cond(Ah) = $(cond(Ah))"
    end

    F = zeros(T, n, n)
    S = zeros(T, n, n)

    metric = one(T)
    it = 0
    while metric > tol && it < maxit
        it += 1
        F = -(Bh + Ch * F) \ Ah
        S = -(Bh + Ah * S) \ Ch

        metric1 = maximum(abs.(Ah + Bh*F + Ch*F*F))
        metric2 = maximum(abs.(Ah*S*S + Bh*S + Ch))
        metric  = max(metric1, metric2)
        # @info "t_iteration step $it, metric = $metric"
    end
    if it == maxit
        @warn "t_iteration did not converge, final metric = $metric"
    end

    eig_F = maximum(abs.(eigvals(F)))
    eig_S = maximum(abs.(eigvals(S)))
    if eig_F > 1/(1-μ) || eig_S >= 1/(1+μ)
        @warn "Conditions of Proposition 3 violated: eig_F=$eig_F, eig_S=$eig_S"
    end

    F_final = F + Iμ                    # Matlab: F = F + I;
    Q_struct = -(B + C * F_final) \ Matrix{T}(I, n, n)  # Matlab: Q = -inv(B + C*F);

    return F_final, Q_struct
    end

    function SolveQZ(Γ0,Γ1,Ψ,Π)
    
    div = 1.0 + 1e-10
    eps = 1e-10
    F = schur!(complex(Γ0),complex(Γ1))
    Lambda, Omega = F.S, F.T
    alpha, beta = F.alpha, F.beta
    Q, Z = adjoint(conj(F.Q)), F.Z

    n = size(Lambda, 1)
    neta = size(Π, 2)

    dLambda = abs.(diag(Lambda))
    dOmega = abs.(diag(Omega))
    dLambda = max.(dLambda,fill(1e-10,size(dLambda))) #to avoid dividing by 0;
    movelast = Bool[(dLambda[i] <= 1e-10) || (dOmega[i] > div * dLambda[i]) for i in 1:n]
    nunstable = sum(movelast)
    nstable = n-nunstable
    iStable = 1:nstable
    iUnstable = (nstable + 1):n

    #Reorder schur to have explosive eigenvalues at the end
    movelastno = fill(false,size(movelast))
    for i in eachindex(movelast)
        movelastno[i] = !movelast[i]
    end
    FS = ordschur!(F, movelastno)
    Lambda, Omega, Q, Z = FS.S, FS.T, FS.Q, FS.Z
    #@show abs.(diag(Lambda))

    gev = hcat(dLambda, dOmega)
    q1 = Q[:,iStable]
    q2 = Q[:,iUnstable]
    q2xΠ = adjoint(q2) * Π
    q2xΨ = adjoint(q2) * Ψ
    q1xΠ = adjoint(q1) * Π
    ndeta1 = min(n - nunstable, neta)
    
    rq2   = rank(q2xΠ)
    rq2q2 = rank([q2xΨ q2xΠ])
    iexist = rq2 == rq2q2
    iunique = rank(Q * Π) == rank(q2xΠ)
    eu = hcat(iexist,iunique)

    #Solve q1xΠ = Phi*q2xΠ by svd decomposition
    #Phi = U1*D1*V1' * V2*inv(D2)*U2
    A2Π = svd(q2xΠ)
    A2Ψ = svd(q2xΨ)
    A1Π = svd(q1xΠ)
    bigevΠ2 = findall(A2Π.S .> eps)
    bigevΨ2 = findall(A2Ψ.S .> eps)
    bigevΠ1 = findall(A1Π.S .> eps)  
    ueta2, deta2, veta2 = A2Π.U[:,bigevΠ2],Matrix(Diagonal(A2Π.S[bigevΠ2])),A2Π.V[:,bigevΠ2]  
    teta, seta, weta = A2Ψ.U[:,bigevΨ2],Matrix(Diagonal(A2Ψ.S[bigevΨ2])),A2Ψ.V[:,bigevΨ2]
    ueta1, deta1, veta1 = A1Π.U[:,bigevΠ1],Matrix(Diagonal(A1Π.S[bigevΠ1])),A1Π.V[:,bigevΠ1]
    Phi =  (ueta1 * deta1 * adjoint(veta1)) * (veta2 * (deta2 \ adjoint(ueta2)))

    #See page 12 of Sims rational expectations document
    L11 = Lambda[iStable,iStable]
    L12 = Lambda[iStable,iUnstable]
    L22 = Lambda[iUnstable,iUnstable]

    O11 = Omega[iStable,iStable]
    O12 = Omega[iStable,iUnstable]
    O22 = Omega[iUnstable,iUnstable]

    Z1 = Z[:,iStable]
    
    #Solve for the effect on lagged variables
    L11inv = LinearAlgebra.inv(L11)
    aux1 = hcat(O11,O12 - Phi*O22) * adjoint(Z)
    aux2 = Z1*LinearAlgebra.inv(L11)
    G1 = real(aux2*aux1)

    #Solve for the effect of exogenous variables (Impact)
    aux3 = vcat(hcat(L11inv, -L11inv*(L12-Phi*L22)),hcat(fill(0.0,(nunstable,nstable)),Matrix(I,nunstable,nunstable)))
    H = Z*aux3
    Impact = real(H * vcat(adjoint(q1) - Phi*adjoint(q2),fill(0.0,(nunstable,size(Ψ,1)))) * Ψ)

    #Solve for the constant 
    #tmat = hcat(Matrix(I,nstable,nstable), -Phi)
    #G0 = vcat(tmat * Lambda, hcat(zeros(nunstable,nstable), Matrix(I,nunstable,nunstable)))
    #G = vcat(tmat * Omega, fill(0.0,(nunstable, n)))
    #G0I = inv(G0)
    #G = G0I * G
    #usix = (nstable + 1):n
    #Ostab = Omega[nstable+1:n,nstable+1:n]
    #Lstab = Lambda[nstable+1:n,nstable+1:n]
    #C = G0I * vcat(tmat * adjoint(Q) * C, (Lstab - Ostab) \ adjoint(q2) * C)

    return eu,G1,Impact
    end

    function TurnABCEtoSims(A,B,C,E)
        HasLead = any((abs.(A) .> 1e-9),dims = 2)
        HasLead = reshape(HasLead,size(A,1))
        Ashift = copy(A)
        Bshift = copy(B)
        Cshift = copy(C)

        Ashift[.!HasLead,:] = B[.!HasLead,:]
        Bshift[.!HasLead,:] = C[.!HasLead,:]
        Cshift[.!HasLead,:] .= 0.0

        #IsLag = findall(any((abs.(Cshift) .> 1e-9),dims=1))
        ##Not sure why I have to use this Linear Indices function, but not using gives me an error in the adjcost case
        IsLag = any((abs.(Cshift) .> 1e-9),dims=1)
        IsLag = LinearIndices(IsLag)[findall(IsLag)]
        n = size(A,1)
        naux = length(IsLag)
        iaux = n+1:n+naux
        G = fill(0.0,(n+naux,n+naux))
        H = fill(0.0,(n+naux,n+naux))

        G[1:n,1:n] = -Ashift
        H[1:n,1:n] = Bshift
        H[1:n,iaux] = Cshift[:,IsLag]

        for (i,col) in enumerate(IsLag)
            G[n+i,n+i] = 1.0
            H[n+i,col] = 1.0
        end

        nEE = length(findall(HasLead))
        EE = fill(0.0,(n+naux,nEE))
        leadeqs = findall(HasLead)
        for (i,lead) in enumerate(leadeqs)
            EE[lead,i] = 1.0
        end
        nE = size(E,2)
        E = vcat(E,fill(0.0,(naux,nE)))
        
        G,H,E,EE = convert(Array{Float64},G),convert(Array{Float64},H),convert(Array{Float64},E),convert(Array{Float64},EE)
        return G,H,E,EE
    end
    using LinearAlgebra

"""
    impresp_sims(G1, Impact, horizon, scale_shocks, scale_y)

Compute impulse responses for a linear RE solution in Sims form.

Model: x_t = G1 * x_{t-1} + Impact * ε_t

Inputs
------
- G1           :: AbstractMatrix{T}    (m_states × m_states)
- Impact       :: AbstractMatrix{T}    (m_states × k_exog)
- horizon      :: Integer              (number of periods)
- scale_shocks :: AbstractVector{T}    (length k_exog, scaling of each shock)
- scale_y      :: Union{T,AbstractVector{T}}
      scalar  → same scaling for all m_states
      vector  → length m_states, per-variable scaling

Output
------
- Resp_mat :: Array{T,3}  (m_states + k_exog, horizon, k_exog)
"""
function impresp_sims(G1::AbstractMatrix{T},
                      Impact::AbstractMatrix{T},
                      horizon::Integer,
                      scale_shocks::AbstractVector{T},
                      scale_y) where {T<:Real}

    m_states, k_exog = size(Impact)
    nz = k_exog

    # 辅助“冲击状态”矩阵（和 matlab 里的 NN 一样）
    NN = zeros(T, nz, nz)

    # 结果：状态 + k_exog 个“冲击自身”变量
    Resp_mat = zeros(T, m_states + k_exog, horizon, k_exog)

    # ---------- 构造 II_lag ----------
    # [ G1   0  ]
    # [ 0   NN  ]
    II_lag = [G1                         zeros(T, m_states, k_exog);
              zeros(T, k_exog, m_states) NN]

    # ---------- 构造 II_contemp ----------
    # aux1 = [ 0(m×m)  Impact ]
    # aux3 = [ 0(k×m)  0(k×k) ]
    aux1 = [zeros(T, m_states, m_states) Impact]
    aux3 = zeros(T, k_exog, m_states + k_exog)
    II_contemp = Matrix{T}(I, m_states + k_exog, m_states + k_exog) + [aux1; aux3]

    # ---------- 变量缩放 ----------
    # scale_y 可以是标量，也可以是长度为 m_states 的向量
    sy_vec = isa(scale_y, Number) ? fill(T(scale_y), m_states) : T.(scale_y)
    @assert length(sy_vec) == m_states

    # 扩展为 (m_states + k_exog)：下面 k_exog 个“冲击自身变量”不缩放
    scale_vec_full = zeros(T, m_states + k_exog)
    scale_vec_full[1:m_states] .= sy_vec

    # ---------- 循环各个冲击 ----------
    for j in 1:k_exog
        Response = zeros(T, m_states + k_exog, horizon)

        # 初始：在“冲击自身”上放 1 单位冲击
        # iVar = m_states + j 对应扩展向量的第 j 个冲击位置
        iVar = m_states + j
        Response[iVar, 1] = one(T)

        # 当期加入 Impact * ε_1
        Response[:, 1] = II_contemp * Response[:, 1]

        # 之后只传播（冲击设为 0）
        M = II_contemp * II_lag    # 提前乘好，减少运算
        for t in 2:horizon
            Response[:, t] = M * Response[:, t-1]
        end

        # shock j 的幅度 × 变量缩放
        s = scale_shocks[j]
        # 这里利用广播: (s .* scale_vec_full) 是长度 m_states+k_exog 的向量
        Response .= (s .* scale_vec_full) .* Response

        Resp_mat[:, :, j] = Response
    end

    return Resp_mat
end
using Plots

"""
    plot_irfs(Resp_mat, n_shocks, vlab; horizonmax=nothing)

在 Julia 中仿照 Matlab 代码绘制 IRFs。

输入
-----
- Resp_mat :: Array{T,3}   (nVarIRF × H × n_shocks)
- n_shocks :: Int          (= n.Shocks)
- vlab     :: Vector{<:AbstractString}  (长度 = nVarIRF)
- horizonmax :: Int (可选)  用于设置 x 轴最大值，默认 = H-1

行为
-----
对每个冲击 ishock 画一张图：
- 前 nVarIRF - n_shocks 行：所有“状态/内生变量”的 IRF
- 额外再画 1 行：该冲击对应的“shock 变量”的 IRF
布局与 Matlab 的 subplot 一致。
"""
function irf_aiy(result; T = 40)
    P     = result.P       # = G1 from SolveQZ
    Q     = result.Q       # = Impact from SolveQZ
    param = result.param
    K_ss  = result.K_ss
    Z_ss  = result.Z_ss

    @unpack nA, nS, α, δ, Lbar = param
    m  = nA * nS            # #assets × #idiosyncratic states
    idxK = 2*m + 1
    idxZ = 2*m + 2

    nx  = size(P, 1)
    IRF = zeros(nx, T)
    IRF[:, 1] = Q[:, 1]     # ε_1 = 1, 之后为 0

    for t in 2:T
        IRF[:, t] = P * IRF[:, t-1]
    end

    # ---- K, Z 的 IRF（百分比）----
    dK = IRF[idxK, :]               # 水平偏差
    dZ = IRF[idxZ, :]

    K_irf = dK ./ K_ss .* 100.0
    Z_irf = dZ ./ Z_ss .* 100.0

    # ---- 在稳态点计算 Y, r, w 及其导数 ----
    Y_ss = Z_ss * K_ss^α * Lbar^(1-α)
    r_ss = α * Y_ss / K_ss - δ
    w_ss = (1-α) * Y_ss / Lbar

    # 导数
    Y_K = α * Z_ss * K_ss^(α-1) * Lbar^(1-α)
    Y_Z = K_ss^α * Lbar^(1-α)

    r_K = α*(α-1) * Z_ss * K_ss^(α-2) * Lbar^(1-α)
    r_Z = α * K_ss^(α-1) * Lbar^(1-α)

    w_K = (1-α)*α * Z_ss * K_ss^(α-1) * Lbar^(-α)
    w_Z = (1-α)    * K_ss^α          * Lbar^(-α)

    # ---- 一阶近似下的 ΔY_t, Δr_t, Δw_t ----
    dY = Y_K .* dK .+ Y_Z .* dZ
    dr = r_K .* dK .+ r_Z .* dZ
    dw = w_K .* dK .+ w_Z .* dZ

    # 百分比 IRF（相对稳态）
    Y_irf = dY ./ Y_ss .* 100.0
    r_irf = dr ./ r_ss .* 100.0
    w_irf = dw ./ w_ss .* 100.0

    # ---- 画图 ----
    tgrid = 1:T
    p1 = plot(tgrid, Z_irf, lw=2, xlabel="t", ylabel="% dev from SS", label="Z")
    p2 = plot(tgrid, K_irf, lw=2, xlabel="t", ylabel="% dev from SS", label="K")
    p3 = plot(tgrid, Y_irf, lw=2, xlabel="t", ylabel="% dev from SS", label="Y")
    p4 = plot(tgrid, r_irf, lw=2, xlabel="t", ylabel="% dev from SS", label="r")
    p5 = plot(tgrid, w_irf, lw=2, xlabel="t", ylabel="% dev from SS", label="w")

    # 自己选布局，我这里给一个 3×2 的例子
    plot(p1, p2, p3, p4, p5, layout=(3,2), size=(900,700))
end
