module MyNumerical_method
export Tauchen, ChebyshevBasis, rouwenhorst
using Parameters,Distributions,Roots, Plots, LinearAlgebra, NLsolve, FastGaussQuadrature, QuadGK, ForwardDiff
    
    function Tauchen(N, ρ, σ; μ = 0.0, K = 3.0)
        σ_x = sqrt(σ^2/(1-ρ^2))
        x_low = μ - K*σ_x
        x_high = μ + K*σ_x
        Δ = (x_high - x_low)/(N-1)
        X = collect(range(x_low,x_high,N))
        Π = fill(0.0, N, N)
        for i in 1:N 
            for j in 2: N-1
                Π[i,j] = cdf(Normal(), (X[j]+Δ/2 - ρ*X[i]-(1-ρ)*μ)/σ) - cdf(Normal(), (X[j]-Δ/2 - ρ*X[i]-(1-ρ)*μ)/σ)
            end
            Π[i,1] = cdf(Normal(),(X[1]+Δ/2 - ρ*X[i]-(1-ρ)*μ)/σ)
            Π[i,N] = 1 - cdf(Normal(), (X[N]-Δ/2-ρ*X[i]-(1-ρ)*μ)/σ)
        end 
        Π = Π./sum(Π, dims = 2)
        return X,Π
    end

    function TauchenHussey(N, ρ, σ; μ = 0.0)
        X, W = gausshermite(N)
        σ_z = σ/(sqrt(1-ρ^2))
        θ = 1/2 + ρ/4
        if ρ >= 0.9
            σ_hat = θ*σ + (1-θ)*σ_z
        else 
            σ_hat = σ
        end
        S = sqrt(2)*σ_hat.*X.+μ
        Π = fill(0.0,N,N)
        τ = σ_hat / σ
        for i in 1:N 
            for j in 1:N
                Π[i,j] = Π[i,j] = exp(X[j]^2- τ^2 * (X[j]-ρ*X[i])^2) *τ *W[j]/sqrt(pi)        
            end
        end
        Π = Π./sum(Π,dims = 2)
        return S, Π
    end

    function rouwenhorst(N, ρ, σ; μ = 0.0)
    q = (ρ+1.0)/2
    nu = ((N-1.0)/(1.0-ρ^2))^(1/2)*σ
    s = collect(range(μ/(1.0-ρ)-nu, μ/(1.0-ρ)+nu, length = N))

    P = [q 1-q; 1-q q]
    for i in 2:N-1
        P = q*[P zeros(i,1); zeros(1,i+1)] + (1-q)*[zeros(i,1) P; zeros(1,i+1)] + (1-q)*[zeros(1,i+1); P zeros(i,1)] + q*[zeros(1,i+1); zeros(i,1) P]
        P[2:i,:]=P[2:i,:]/2
    end
    return s, P
    end

    function ChebyshevNodes(N_nodes,a,b)
        X = -cos.((2.0 .*collect(1:N_nodes).-1.0).*pi ./(2.0*N_nodes) )
        Y = (X.+1.0 ).*(b-a)./2.0 .+a
        return X,Y
    end

    function ChebyshevBasis(N_order, Z)
    # This is code for ChebyshevBasis
    # N_order : The order of Polynominal
    # Z : The location of nodes
    id = 1
    T = fill(0.0, 1,N_order)
        while id < N_order + 1
            θ = acos(Z)
            T[id] = cos((id - 1) * θ)
            id = id + 1
        end   
        return T     
    end

    function Rootfinding_Newton(f,x_0)
        #这是一个基于牛顿法计算零根的函数
        #f 满足f(x) = 0
        # x_0 是初始猜测
        J(x) = ForwardDiff.jacobian(f,x)
        ϵ = 1E-12
        weight = 0.2
        iter = 1
        maxit = 10000
        crit = 1
        x = x_0
        while crit>ϵ && iter<maxit
            x_new = x - J(x)\f(x)
            crit = maximum(abs.(x_new .- x))
            iter = iter + 1
            x .= weight * x_new + (1 - weight)*x
        end
        return x
    end

    function Rootfinding_Bisection(f,a,b; ϵ = 1e-8, maxit = 500)
        s = sign(f(a))
        x_midpoint = (a+b)/2
        step = (b-a)/2
        #ϵ = 1E-8
        iter = 0
        while step>ϵ && iter < maxit
            step = step/2
            if s == sign(f(x_midpoint))
                x_midpoint = x_midpoint + step
            else
                x_midpoint = x_midpoint - step
            end
            iter = iter +1
        end
        println("Tol is $step,","Iteration is $iter")
        return x_midpoint
    end

    function Rootfinding_NewtonRaphson(f, x_0; ϵ = 1E-12, maxit = 10000, weight = 0.2)
        # 牛顿-拉夫森算法求解无约束最小化问题，注意推导使用了一阶条件！
        J(x) = ForwardDiff.jacobian(f,x)
        H(x) = ForwardDiff.hessian(f,x)
        crit = 1
        iter = 1
        x = x_0
        while iter<maxit && crit>ϵ
            x_new = x - H(x)\J(x)
            crit = maximum(abs.(x_new .- x))
            iter = iter + 1
            x .= weight*x_new + (1 - weight )*x
        end
        return x
    end

    function linearGridIntp!(xgrid::Array{Float64,1}, ygrid::Array{Float64,1}, xtarget::Array{Float64,1}, ytarget::Array{Float64,1})

	# check if is sorted
	@assert (issorted(xgrid) & issorted(xtarget))

	xg = length(xgrid); #xt = length(xtarget)
	xi  = 1; xlow = xgrid[1]; xhigh = xgrid[2] # initialize
	@inbounds for (it, xtval) in enumerate(xtarget)
                    while xi < xg - 1 # find grid point
                        if xhigh >= xtval; break; end
                        xi += 1
                        xlow = xhigh
                        xhigh = xgrid[xi + 1]
                    end
                    xprob = (xhigh - xtval) / (xhigh - xlow)
                    ytarget[it] = ygrid[xi] * xprob + (1.0 - xprob) * ygrid[xi+1]
	            end
	return ytarget
    end

    function coordGridIntp!(xgrid::AbstractVector, xtarget::AbstractVector, 		ibelow::AbstractVector{Int}, iweight::AbstractVector; robust = false)
            xg = length(xgrid); #xt = length(xtarget)
            xi  = 1; xlow = xgrid[1]; xhigh = xgrid[2] # initialize
            @inbounds for (it, xtval) in enumerate(xtarget)
                while xi < xg - 1 # find grid point
                    if xhigh >= xtval; break; end
                    xi += 1
                    xlow = xhigh
                    xhigh = xgrid[xi + 1]
                    end
                iweight[it] = (xhigh - xtval) / (xhigh - xlow) #  careful with last point. if xhigh<xtval might run into problems.
                ibelow[it] = xi
            end
            if robust == true; iweight =  min.(max.(iweight,0.0), 1.0); end
        return ibelow,  iweight
    end

    function linearGridIntp!(xgrid::Array{Float64,1}, ygrid::Array{Float64,1}, xtarget::Array{Float64,1}, ytarget::Array{Float64,1})

	# check if is sorted
	@assert (issorted(xgrid) & issorted(xtarget))

	xg = length(xgrid); #xt = length(xtarget)
	xi  = 1; xlow = xgrid[1]; xhigh = xgrid[2] # initialize
	@inbounds for (it, xtval) in enumerate(xtarget)
		while xi < xg - 1 # find grid point
			if xhigh >= xtval; break; end
			xi += 1
			xlow = xhigh
            xhigh = xgrid[xi + 1]
		end
		xprob = (xhigh - xtval) / (xhigh - xlow)
		ytarget[it] = ygrid[xi] * xprob + (1.0 - xprob) * ygrid[xi+1]
	end
	return ytarget
    end
    """
    通用的 coordGridIntp!，支持 ForwardDiff.Dual 等任意 Number。
    给定单调递增网格 Grid 和查询点向量 X，返回:
    - id_L[i] : 使得 Grid[id_L[i]] ≤ X[i] ≤ Grid[id_L[i]+1]
    - weight_L[i] : 落在左端点上的权重 (右端点权重 = 1 - weight_L[i])
    """
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


    # function that calculates gini (see wikipedia formula):
	function calGini2(dist, values)
		cumulative = cumsum(dist.*values)/sum(dist.*values) # cum distribution in %
		B = sum(cumulative.*dist) # total area below lorenz curve
		# with cns distribution this should be 0.5, but since it is a histogram it may deviate a bit
		AreaBelow45 = sum(dist.*(cumsum(dist)./sum(dist))) # A + B
		gini = (AreaBelow45 - B)/AreaBelow45
		return gini
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
    function lorenz_curve(dist, value; doPlot=true)
    # value: 资产/收入格点 a_j（允许未排序）
    # dist : 对应权重 m_j（概率或计数，允许未归一化）
    # doPlot: 是否画图（默认 true）
    
    value = vec(value)  # 确保是向量
    dist = vec(dist)
    
    # 数据验证
    @assert all(isfinite, value) && all(isfinite, dist) "NaN/Inf detected"
    @assert all(dist .>= 0) "dist 必须非负"
    
    # 排序并同步权重
    sorted_indices = sortperm(value)
    value = value[sorted_indices]
    dist = dist[sorted_indices]
    
    M = sum(dist)        # 总人口
    W = sum(dist .* value) # 总财富
    @assert M > 0 && W >= 0 "权重或总财富异常"
    
    # 累计人口份额与累计财富份额
    p = cumsum(dist) / M
    L = cumsum(dist .* value) / W
    
    # 为了从原点出发，拼上 (0,0)
    p = [0; p]
    L = [0; L]
    
    # 基尼系数 G = 1 - 2*∫_0^1 L(p) dp （梯形积分）
    area = trapz(p, L)
    G = 1 - 2 * area
    
    if doPlot
       p1 = plot(p, L, linewidth=2, label="Lorenz Curve")
        plot!([0, 1], [0, 1], linestyle=:dash, color=:black, linewidth=1, label="Equality Line")
        xlabel!("Cumulative population share")
        ylabel!("Cumulative wealth share")
        title!("Lorenz Curve (Gini = $(round(G, digits=3)))")
        xlims!(0, 1)
        ylims!(0, 1)
        plot(p1)
        #grid!(true)
        #legend()
    end
    
    #return (p, L, G)
    end

    # 梯形积分辅助函数
    function trapz(x, y)
        n = length(x)
        @assert length(y) == n "x和y长度必须相同"
        
        integral = 0.0
        for i in 2:n
            dx = x[i] - x[i-1]
            integral += dx * (y[i] + y[i-1]) / 2
        end
        return integral
    end
    
end
