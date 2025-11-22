module MyNumerical_method
export Tauchen, ChebyshevBasis
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

    function ChebyshevNodes(N_nodes)
        X = -cos.((2 .*collect(1:N_nodes).-1).*pi ./(2*N_nodes) )
        return X
    end

    function ChebyshevBasis(N_order, Z)
    # This is code for ChebyshevBasis
    # N_order : The order of Polynominal
    # Z : The location of nodes
    id = 1
    T = fill(0.0, N_order)
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

    function Rootfinding_Bisection(f,a,b)
        s = sign(f(a))
        x_midpoint = (a+b)/2
        step = (b-a)/2
        ϵ = 1E-12
        while step>ϵ
            step = step/2
            if s == sign(f(x_midpoint))
                x_midpoint = x_midpoint + step
            else
                x_midpoint = x_midpoint - step
            end
        end
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

function coordGridIntp!(xgrid::Array{Float64,1}, xtarget::Array{Float64,1}, 		ibelow::Array{Int64,1}, iweight::Array{Float64,1}; robust = false)
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
	return iweight, ibelow
end

end
