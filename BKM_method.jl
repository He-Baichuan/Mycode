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
    τ_l0 = 0.0, #第0期的劳动税率
    τ_l1 = 0.0, #第一期的劳动税率
    τ_k0 = 0.0, #第零期的资本税率
    τ_k1 = 0.0, #第一期的资本税率
    G0 = 0.0, #第0期的政府支出比率
    G1 = 0.0, #第一期的政府支出比率
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
    Zbar = exp(Zbar)
    Shock_Z[1] = log(Zbar); Shock_Z[T] = log(Zbar); Shock_Z[2] = ShockSize
    for t in 3:T-1
        Shock_Z[t] = ρ_z * Shock_Z[t-1]
    end
    Shock_Z = exp.(Shock_Z)
    #Zbar = exp(Zbar)

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
    """
    利用内生格点法求解家庭消费储蓄问题
    """
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
        
            for i in 1:nA, j in 1:nS 
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
    """
    young(2009)+MC 求解不变分布
    """
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
    #总冲击影响总量和价格，并继续影响家庭决策
    r = Zbar*α *Kguess^(α-1)*Lbar^(1-α)-δ
    w = Zbar*(1-α)*Kguess^α*Lbar^(-α)
    Y = Zbar*Kguess^α*Lbar^(1-α)
    G = Y*Gexp
    transfer = τl*w*Lbar + τk*r*Kguess - G #政府预算平衡
    decision = SolveHHprobEGM(param, r, w, τl, τk, transfer)
    dsn = SolveInvariantDist(param, decision)
    gAA = repeat(gridA,1,nS)
    Ea = sum(dsn.*gAA)
    #excess = 2*(Ea - Kguess)/(Ea + Kguess)
    excess = Ea - Kguess
    return (excess, dsn, decision, w, r, transfer, Ea) 
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
    return (decision, dsn, w, r, transfer, eqm_K, Ea)
end

function BackwardInductionSeq(param, rseq, wseq, TransSeq, PolSSend)
    """
    PolSSend : 新稳态的政策函数（也是逆向归纳法的初始条件）
    同样利用内生格点法从后往前计算每一期的政策函数
    """
    @unpack nA, nS, gridA, gridS, transS, β, γ, δ, T, τk_seq, τl_seq = param 
    gAA = repeat(gridA,1,nS)
    gSS = repeat(gridS',nA,1)

    cpseq = zeros(nA,nS,T)
    apseq = zeros(nA,nS,T)

    apseq[:,:,T] = PolSSend.ap 
    cpseq[:,:,T] = PolSSend.cp 

    for t in T-1:-1:2
        RHS = β*(1.0+rseq[t+1]*(1-τk_seq[t+1]))*cpseq[:,:,t+1].^(-γ)*transS'
        c = RHS.^(-1.0/γ)
        A_endo = (c + gAA - wseq[t].*gSS*(1-τl_seq[t]).-TransSeq[t])/(1.0 + rseq[t]*(1.0 - τk_seq[t]))
        ap_new = zeros(nA,nS)
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
        apseq[:,:,t] = ap_new
        cpseq[:,:,t] = (1.0+rseq[t]*(1.0 - τk_seq[t])).*gAA .+ wseq[t]*(1.0-τl_seq[t]).*gSS .+ TransSeq[t] - ap_new
    end
    return (apseq = apseq, cpseq = cpseq)
end

function SolveDistSeq(param, apseq, dsnbeg)
    @unpack nA, nS, gridA, transS, T = param 
    """
    基于政策函数序列从前向后计算分布序列
    """
    dsnseq = zeros(nA,nS,T)
    dsnseq[:,:,1] = dsnbeg # 注意到分布是期末确定的，因此前两期的分布是给定的
    dsnseq[:,:,2] = dsnbeg

    for t in 2:T-2 
        id_L = fill(0, nA, nS)
        weight_L = zeros(nA, nS)
        for j in 1:nS 
            id_L[:,j], weight_L[:,j] = MyNumerical_method.coordGridIntp!(gridA, apseq[:,j,t],id_L[:,j], weight_L[:,j],robust = true)
        end
        dsn_new = zeros(nA,nS)
        for i in 1:nA, j in 1:nS 
            if dsnseq[i,j,t] >0.0
               dsn_new[id_L[i,j],j] += weight_L[i,j]*dsnseq[i,j,t]
               dsn_new[id_L[i,j]+1,j] += (1.0 - weight_L[i,j])*dsnseq[i,j,t] 
            end
        end
        dsnseq[:,:,t+1] = dsn_new*transS
    end
    return dsnseq
end

function ExcessDemandTrans(param, solSS0, solSS1, Kseq)
    @unpack gridA, nS, Lbar, α, δ, T, Shock_Z, Gseq, τl_seq, τk_seq = param 

    rseq = zeros(T)
    wseq = zeros(T)
    rseq[1] = solSS0[4]
    rseq[T] = solSS1[4]
    wseq[1] = solSS0[3]
    wseq[T] = solSS1[3]
    for t in 2:T-1 
        rseq[t] = Shock_Z[t]*α*(Kseq[t]/Lbar)^(α - 1)-δ
        wseq[t] = Shock_Z[t]*(1-α)*(Kseq[t]/Lbar)^(α)
    end

    Y = Shock_Z.*Kseq.^α.*Lbar^(1-α)
    G = Y.*Gseq
    TransSeq = rseq.*Kseq.*τk_seq .+ wseq.*τl_seq.*Lbar - G 

    apseq,cpseq = BackwardInductionSeq(param, rseq, wseq, TransSeq, solSS1[1])
    apseq[:,:,1] = solSS0[1].ap
    cpseq[:,:,1] = solSS0[1].cp

    dsnseq = SolveDistSeq(param, apseq, solSS0[2])
    dsnseq[:,:,T] = solSS1[2]

    Kseq_next = zeros(T)
    Kseq_next[1] = solSS0[6]
    Kseq_next[2] = solSS0[6]
    Kseq_next[T] = solSS1[6]

    for t in 3:T-1
        Kseq_next[t] = sum(dsnseq[:,:,t].*repeat(gridA,1,nS))
    end
    return Kseq_next, dsnseq, apseq, cpseq, wseq, rseq, TransSeq
end

function SolveModelTrans(param)
    T = param.T 
    VaratTime = [param.Shock_Z[1], param.τl_seq[1], param.τk_seq[1], param.Gseq[1]]
    solSS0 = SolveModel(param, VaratTime)
    VaratTime = [param.Shock_Z[T], param.τl_seq[T], param.τk_seq[T], param.Gseq[T]]
    solSS1 = SolveModel(param, VaratTime)
    #初始猜测资本总量是线性增长的
    Kseq = zeros(T)
    Kseq[1] = solSS0[6]
    Kseq[T] = solSS1[6]
    d = (Kseq[T]-Kseq[1])/(T-1)
    for i in 2:T-1
        Kseq[i] = Kseq[1]+d*(i-1)
    end
    maxit = 5000
    tol = 1e-7
    damp = 0.2

    for iter in 1:maxit
        Kseq_next, = ExcessDemandTrans(param, solSS0, solSS1, Kseq)

        crit = maximum(abs.(Kseq_next - Kseq))
        Kseq = damp*Kseq_next + (1.0 - damp)*Kseq
        if crit<tol
            println("done!")
            break
        end
        if iter == maxit
            println("No done")
        end
    end
    res = ExcessDemandTrans(param, solSS0, solSS1, Kseq)
    return res
end
res = SolveModelTrans(SettingPar_TransitionDynamics(;T = 180))
param = SettingPar_TransitionDynamics(;T = 180)
@unpack Shock_Z, T, gridA, nS = param 
gridT = collect(range(1,T,T))
gini = zeros(T,3)
Share_const = zeros(T)

for t in 1:T 
    dsn = res[2][:,:,t]
    asset_dist = sum(dsn, dims = 2)[:]
    gini[t,1] = MyNumerical_method.calGini2(asset_dist,gridA)
    gini[t,2] = MyNumerical_method.calGini2(dsn[:,1]./sum(dsn[:,1]),gridA)
    gini[t,3] = MyNumerical_method.calGini2(dsn[:,end]./sum(dsn[:,end]),gridA)
    Share_const[t] = dsn[1]
end

Zper = log.(Shock_Z)*100
Kper = (log.(res[1]).-log(res[1][1])).*100
rpp = (res[6].- res[6][1]).*100
wper = (log.(res[5]).-log(res[5][1])).*100
share_pp = (Share_const.-Share_const[1]).*100
gini1_pp =  ( gini[:,1]  .-  gini[1,1])*100
gini2_pp =  ( gini[:,2]  .-  gini[1,2])*100
gini3_pp =  ( gini[:,3]  .-  gini[1,3])*100
p1 = plot(gridT, Zper, label="Z", lw = 3, ylabel = "%", xlabel = "t", xlims = (0, 180), ylims = (0,1))
p2 = plot(gridT, Kper, label="K", lw = 3, ylabel = "%", xlabel = "t", xlims = (0, 180), ylims = (0,1))
p3 = plot(gridT, rpp, label="r", lw = 3, ylabel = "p.p", xlabel = "t", xlims = (0, 180))
p4 = plot(gridT, wper, label="w", lw = 3, ylabel = "%", xlabel = "t", xlims = (0, 180), ylims = (0,1))
pp1 = plot(p1, p2, p3, p4, layout = (2, 2), size = (900,800))
p5 = plot(gridT[1:end-1], gini1_pp[1:end-1], title="Gini Wealth (All)", lw = 3, ylabel = "p.p", xlabel = "t", xlims = (0,180))
p6 = plot(gridT[1:end-1], gini2_pp[1:end-1], title="Gini Wealth (Smin)", lw = 3, ylabel = "p.p", xlabel = "t", xlims = (0,180))
p7 = plot(gridT[1:end-1], gini3_pp[1:end-1], title="Gini Wealth (Smax)", lw = 3, ylabel = "p.p", xlabel = "t", xlims = (0,180))
p8 = plot(gridT[1:end-1], share_pp[1:end-1], title="Share Constrained", lw = 3, ylabel = "p.p", xlabel = "t", xlims = (0,180))
pp2 = plot(p5, p6, p7, p8, layout = (2, 2), legends = false, size = (900,800))
using Random
using Distributions
@unpack ρ_z, Zbar = param 
simul = 500
Random.seed!(20251122)

K_irf = res[1][2:T]
Kss = res[1][1]
n_shocks = length(K_irf)


shocks = rand(Normal(0,1), simul)
log_z_simul = zeros(length(shocks))
log_z_simul[1] = ρ_z*log(Zbar) + shocks[1]
for t in 2:length(shocks)
    log_z_simul[t] = ρ_z*log_z_simul[t-1]+shocks[t]
end
deriv = K_irf./Kss .- 1.0 
function BPM_path(deriv, shocks, simul, n_shocks)
    K_simul = zeros(simul)
    for t in 1:simul
    tmp = 0.0
    Jt = min(t, n_shocks)      # 最多只能用到 t 期以前的那些 IRF
        for j in 1:Jt
            # j 对应 ψ_{j-1}，冲击是 ε_{t-j+1}
            tmp += deriv[j] * shocks[t - j + 1]
        end
    
        K_simul[t] = 100.0 * tmp
    end
    return K_simul
end
    K_simul = BPM_path(deriv, shocks, simul, n_shocks)

	gSimul = collect(range(1,simul,step=1))
	
	pK = plot(gSimul, K_simul, lw = 3, ylabel = "% deviations of SS", xlabel = "t", label = "K")
	pZ = plot(gSimul, log_z_simul, lw = 3, ylabel = "% deviations of SS", xlabel = "t", label = "log Z")
	plot(pK, pZ, layout = (2, 1))
    
param2 = SettingPar_TransitionDynamics(ShockSize = 0.0, τ_l1 = 0.2, T = 160)
res2 = SolveModelTrans(param2)

	τlseq = param2.τl_seq*100 
	Kper2 = @. (log(res2[1]) - log(res2[1][1]))*100
	rpp2 = @. (res2[6] - res2[6][1])*100
	wper2 = @. (log(res2[5]) - log(res2[5][1]))*100
    T = param2.T
	gT = collect(range(1,T))
	p21 = plot(gT, τlseq, title="Labor Tax", lw = 3, ylabel = "%", xlabel = "t")
	p22 = plot(gT, Kper2, title="K", lw = 3, ylabel = "%", xlabel = "t")
	p23 = plot(gT, rpp2, title="r", lw = 3, ylabel = "p.p", xlabel = "t")
	p24 = plot(gT, wper2, title="w", lw = 3, ylabel = "%", xlabel = "t")
	plot(p21, p22, p23, p24, layout = (2, 2), legends = false)
