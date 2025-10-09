include("Numerical_method.jl")
#parameter setting
N = 3
σ = 4.0
A = [2.0; 1.0; 1.0]
L = [1.0; 2.0; 4.0]
τ = [1.0 2.0 2.0; 2.0 1.0 2.0; 2.0 2.0 1.0]
T = fill(0.0, N, N)

function eqm_update(w,x,T)
    κ = τ .* (1 .+ T)
    K = fill(0.0, N, N)
    # build K elementwise
    for i in 1:N, j in 1:N
        K[i,j] = ( (w[i]/A[i]) * κ[i,j] )^(1-σ)
    end

    P = fill(0.0, N)
    # sum over origin i for each destination n (column n)
    for n in 1:N
        P[n] = sum(K[:, n])^(1/(1-σ))
    end

    welfare = x ./ P

    λ = fill(0.0, N, N)
    for i in 1:N, j in 1:N
        λ[i,j] = K[i,j] * P[j]^(σ - 1)
    end

    x_rep = [x'; x'; x']
    X_in = λ .* x_rep

    w_new = sum(X_in ./ (1 .+ T); dims =2) ./ L
    x_new = w .* L + sum(X_in .* (T ./ (1 .+ T)); dims = 2)

    wnum = sum(w_new)
    w_new = w_new ./ wnum
    x_new = x_new ./ wnum
    P = P ./ wnum

    return w_new, x_new, P, welfare, λ
end

function eqm_iter(T; ϵ = 1E-6, maxit = 10000, weight = 0.2)
    crit = 1.0
    iter = 1
    w = fill(1.0 / N, N)
    x = copy(w)
    P = fill(0.0,N)
    welfare = copy(P)
    λ = fill(0.0,N,N)
    while crit > ϵ && iter < maxit
        w_new, x_new, P, welfare, λ = eqm_update(w, x, T)
        r = [w; x]
        r_new = [w_new; x_new]
        crit = maximum(abs.(r_new - r))
        iter += 1
        w .= weight * w_new + (1 - weight) * w
        x .= weight * x_new + (1 - weight) * x
    end
    return w, x, P, welfare, λ
end

sol = eqm_iter(T)
T_smul = copy(T)
T_smul[2,1] = 0.5
T_smul[3,1] = 0.5
sol2 = eqm_iter(T_smul)
println(sol[4])
println(sol2[4])

