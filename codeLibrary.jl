using LinearAlgebra


function LinearTimeIiteration(A,B,C,D; maxit = 10000, tol = 1e-8)
    P0 = fill(0.0, size(A))
    S0 = fill(0.0, size(C))
    for it in 1:maxit
        P1 = -(A*P0+B)\C 
        S1 = -(C*S0+B)\A 
        if maximum(abs.(P1 - P0))< tol
            break
        end
        P0 = P1
        S0 = S1
    end
    Q = -(A*P0+B)\D 
    eigen_P = eigen(P0)
    eigen_S = eigen(S0)
    if maximum(abs.(eigen_P.values)) > 1.0
        error("Non existence")
    end
    if maximum(abs.(eigen_S.values)) > 1.0
        error("No stable equilibrium")
    end
    return P0, Q
end

function interpEGM(
    xgrid::AbstractVector,
    ygrid::AbstractVector,
    x::Real
)

    n = length(xgrid)

    j = clamp(
        searchsortedlast(xgrid, x),
        1,
        n - 1
    )

    x1 = xgrid[j]
    x2 = xgrid[j+1]

    y1 = ygrid[j]
    y2 = ygrid[j+1]

    ω = (x - x1) / (x2 - x1)

    y = y1 + ω * (y2 - y1)

    return y, j
end