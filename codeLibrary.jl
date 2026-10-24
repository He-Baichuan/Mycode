using LinearAlgebra


function LinearTimeIteration(A,B,C,D; maxit = 10000, tol = 1e-8)
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

# ============================================================
# Numerical rank information
# ============================================================

function numerical_rank_info(
    M;
    rtol = 1e-10,
    atol = 0.0
)

    s = svdvals(Matrix(M))

    if isempty(s)
        return (
            rank = 0,
            fullrank = false,
            σmax = 0.0,
            σmin = 0.0,
            cutoff = 0.0
        )
    end

    σmax = maximum(s)

    cutoff =
        max(
            atol,
            rtol * σmax
        )

    r =
        count(x -> x > cutoff, s)

    σmin =
        minimum(s)

    return (
        rank = r,
        fullrank = (r == min(size(M)...)),
        σmax = σmax,
        σmin = σmin,
        cutoff = cutoff
    )
end

# ============================================================
# Rank-revealing QR
#
# A ≈ U * L
#
# U : n × r
# L : r × n
# ============================================================

function lowrank_qr(
    A;
    rtol = 1e-10,
    atol = 0.0
)

    A0 = Matrix(A)

    n,m =
        size(A0)

    F =
        qr(
            A0,
            ColumnNorm()
        )

    R =
        Matrix(F.R)

    Q =
        Matrix(F.Q)

    p =
        F.p


    d =
        abs.(diag(R))


    if isempty(d) ||
       maximum(d) == 0

        T =
            eltype(A0)

        return (
            U = zeros(T,n,0),
            L = zeros(T,0,m),
            rank = 0,
            error = 0.0
        )
    end


    cutoff =
        max(
            atol,
            rtol * maximum(d)
        )

    r =
        count(
            x -> x > cutoff,
            d
        )


    U =
        Q[:,1:r]

    Lp =
        R[1:r,:]


    # A[:,p] = Q*R
    # undo column permutation
    L =
        zeros(
            eltype(A0),
            r,
            m
        )

    L[:,p] =
        Lp


    scale =
        max(
            opnorm(A0,Inf),
            eps(Float64)
        )


    qr_error =
        opnorm(
            A0-U*L,
            Inf
        ) / scale


    return (
        U = U,
        L = L,
        rank = r,
        error = qr_error
    )
end
function _LITx_anchored(
    A,
    B,
    C,
    D,
    U,
    L,
    P0;
    maxit = 10_000,
    tol = 1e-10,
    damping = 1.0,
    verbose = true
)

    n =
        size(A,1)

    r =
        size(U,2)


    # --------------------------------------------------------
    # Anchor
    #
    # W0 = L*P0
    #
    # M0 = B + A*P0
    #    = B + U*W0
    # --------------------------------------------------------

    W0 =
        L * P0

    M0 =
        B +
        U * W0


    FM0 =
        lu(Matrix(M0))


    # --------------------------------------------------------
    # Anchored Woodbury objects
    #
    # G0 = -M0^{-1} C
    #
    # H0 =  M0^{-1} U
    # --------------------------------------------------------

    G0 =
        -(FM0 \ Matrix(C))

    H0 =
        FM0 \ U


    LG0 =
        L * G0

    LH0 =
        L * H0


    T =
        promote_type(
            eltype(A),
            eltype(B),
            eltype(C),
            eltype(D)
        )


    Ir =
        Matrix{T}(I,r,r)


    # start exactly from the anchor
    W =
        copy(W0)


    err =
        Inf

    nit =
        0

    converged =
        false


    for it in 1:maxit

        # ----------------------------------------------------
        # ΔW = W - W0
        # ----------------------------------------------------

        ΔW =
            W - W0


        # ----------------------------------------------------
        # Small r × r system
        #
        # Msmall = I + ΔW H0
        # ----------------------------------------------------

        Msmall =
            Ir +
            ΔW*H0


        Y =
            Msmall \ (ΔW*G0)


        # ----------------------------------------------------
        # Exact anchored LIT update
        #
        # Wnew =
        #
        # L G0
        # - L H0 (I+ΔW H0)^(-1) ΔW G0
        # ----------------------------------------------------

        Wcandidate =
            LG0 -
            LH0*Y


        Wnew =
            damping .* Wcandidate +
            (1.0-damping) .* W


        err =
            maximum(
                abs.(Wnew-W)
            ) /
            max(
                1.0,
                maximum(abs.(Wnew))
            )


        W =
            Wnew

        nit =
            it


        if verbose &&
           (it == 1 || it % 10 == 0)

            println(
                "LITx iteration ",
                it,
                ": error = ",
                err
            )
        end


        if err < tol

            converged =
                true

            break
        end
    end


    if !converged
        error(
            "Anchored LITx did not converge after $maxit iterations."
        )
    end


    # ========================================================
    # Recover full P
    # ========================================================

    ΔW =
        W - W0

    Msmall =
        Ir +
        ΔW*H0


    P =
        G0 -
        H0 * (
            Msmall \
            (ΔW*G0)
        )


    # ========================================================
    # Shock impact
    #
    # Q = -(A P+B)^(-1)D
    # ========================================================

    J0 =
        -(FM0 \ Matrix(D))


    Qshock =
        J0 -
        H0 * (
            Msmall \
            (ΔW*J0)
        )


    # ========================================================
    # Diagnostics
    # ========================================================

    RP =
        A*P*P +
        B*P +
        C

    RQ =
        (A*P+B)*Qshock +
        D


    return (
        P = P,
        Q = Qshock,
        iterations = nit,
        error = err,
        residual_P = maximum(abs.(RP)),
        residual_Q = maximum(abs.(RQ))
    )
end


# ============================================================
# Generic LITx
#
# 1. QR compress A
# 2. Check B
#
# If B full rank:
#     ordinary QR-Woodbury LITx
#
# If B singular:
#     regularize B
#     solve regularized system -> Preg
#     construct original M0 = B + A*Preg
#
#     if M0 full rank:
#         anchored LITx on ORIGINAL model
#
#     else:
#         fallback to ordinary LIT
# ============================================================

function LinearTimeIterationX(
    A,
    B,
    C,
    D;
    maxit = 10_000,
    tol = 1e-10,
    ranktol = 1e-10,
    rankatol = 0.0,
    damping = 1.0,

    # regularization
    reg_lambda = 1e-8,
    reg_growth = 10.0,
    reg_max_tries = 5,
    reg_tol = 1e-8,

    verbose = true
)

    n =
        size(A,1)


    @assert size(A) == (n,n)
    @assert size(B) == (n,n)
    @assert size(C) == (n,n)
    @assert size(D,1) == n


    # ========================================================
    # 1. QR compression of future Jacobian A
    # ========================================================

    qrA =
        lowrank_qr(
            A;
            rtol = ranktol,
            atol = rankatol
        )


    U =
        qrA.U

    L =
        qrA.L

    rA =
        qrA.rank


    if verbose

        println(
            "LITx numerical rank(A) = ",
            rA,
            " / ",
            n
        )

        println(
            "QR approximation error = ",
            qrA.error
        )
    end


    # ========================================================
    # Special case A ≈ 0
    # ========================================================

    if rA == 0

        FB =
            lu(Matrix(B))

        P =
            -(FB \ Matrix(C))

        Q =
            -(FB \ Matrix(D))


        return P,Q,(
            method = :static,
            rankA = 0,
            regularized = false
        )
    end


    # ========================================================
    # 2. Check B rank
    # ========================================================

    infoB =
        numerical_rank_info(
            B;
            rtol = ranktol,
            atol = rankatol
        )


    if verbose

        println(
            "rank(B) = ",
            infoB.rank,
            " / ",
            n
        )

        println(
            "σmin(B) = ",
            infoB.σmin
        )
    end


    # ========================================================
    # CASE 1:
    # B is full rank
    #
    # ordinary QR-Woodbury LITx
    # ========================================================

    if infoB.fullrank

        if verbose
            println(
                "B is full rank -> standard QR-Woodbury LITx."
            )
        end


        P0 =
            zeros(
                promote_type(
                    eltype(A),
                    eltype(B),
                    eltype(C)
                ),
                n,n
            )


        sol =
            _LITx_anchored(
                A,
                B,
                C,
                D,
                U,
                L,
                P0;
                maxit = maxit,
                tol = tol,
                damping = damping,
                verbose = verbose
            )


        info = (
            method = :litx,
            rankA = rA,
            qr_error = qrA.error,
            rankB = infoB.rank,
            regularized = false,
            M0_fullrank = true,
            iterations = sol.iterations,
            error = sol.error,
            residual_P = sol.residual_P,
            residual_Q = sol.residual_Q
        )


        return sol.P,sol.Q,info
    end


    # ========================================================
    # CASE 2:
    # B is rank deficient
    #
    # use regularization ONLY to generate an anchor Preg
    # ========================================================

    if verbose

        println(
            "B is rank deficient."
        )

        println(
            "Start regularized initialization..."
        )
    end


    scaleB =
        max(
            opnorm(B,Inf),
            1.0
        )


    Preg =
        nothing

    λused =
        NaN

    reg_success =
        false


    λ =
        reg_lambda


    for k in 1:reg_max_tries

        Breg =
            Matrix(B) +
            λ * scaleB * I


        infoBreg =
            numerical_rank_info(
                Breg;
                rtol = ranktol,
                atol = rankatol
            )


        if verbose

            println(
                "regularization try ",
                k,
                ": λ = ",
                λ,
                ", rank(Breg) = ",
                infoBreg.rank,
                "/",
                n
            )
        end


        if !infoBreg.fullrank

            λ *= reg_growth

            continue
        end


        # ----------------------------------------------------
        # Solve REGULARIZED model only to obtain Preg
        # ----------------------------------------------------

        Pzero =
            zeros(
                promote_type(
                    eltype(A),
                    eltype(B),
                    eltype(C)
                ),
                n,n
            )


        try

            regsol =
                _LITx_anchored(
                    A,
                    Breg,
                    C,
                    D,
                    U,
                    L,
                    Pzero;
                    maxit = maxit,
                    tol = reg_tol,
                    damping = damping,
                    verbose = false
                )


            Preg =
                regsol.P

            λused =
                λ

            reg_success =
                true


            if verbose

                println(
                    "Regularized initialization succeeded with λ = ",
                    λused
                )
            end


            break

        catch err

            if verbose

                println(
                    "Regularized LITx failed at λ = ",
                    λ,
                    ". Increasing λ."
                )
            end

            λ *= reg_growth
        end
    end


    # ========================================================
    # Regularization itself failed
    #
    # fallback directly to LIT
    # ========================================================

    if !reg_success

        if verbose

            println(
                "Regularized initialization failed."
            )

            println(
                "Fallback -> ordinary LIT."
            )
        end


        P,Q =
            LinearTimeIteration(
                A,
                B,
                C,
                D;
                maxit = maxit,
                tol = tol
            )


        info = (
            method = :lit_fallback,
            rankA = rA,
            qr_error = qrA.error,
            rankB = infoB.rank,
            regularized = true,
            regularization_success = false,
            M0_fullrank = false
        )


        return P,Q,info
    end


    # ========================================================
    # 3. Return to ORIGINAL model
    #
    # M0 = B + A*Preg
    # ========================================================

    M0 =
        B +
        A*Preg


    infoM0 =
        numerical_rank_info(
            M0;
            rtol = ranktol,
            atol = rankatol
        )


    if verbose

        println(
            "rank(M0 = B + A*Preg) = ",
            infoM0.rank,
            " / ",
            n
        )

        println(
            "σmin(M0) = ",
            infoM0.σmin
        )

        println(
            "M0 full rank = ",
            infoM0.fullrank
        )
    end


    # ========================================================
    # CASE 2A:
    # M0 is full rank
    #
    # Use Preg only as anchor,
    # then solve ORIGINAL system exactly.
    # ========================================================

    if infoM0.fullrank

        if verbose

            println(
                "M0 is full rank."
            )

            println(
                "Switching back to ORIGINAL model "
                * "with anchored QR-Woodbury LITx."
            )
        end


        sol =
            _LITx_anchored(
                A,
                B,
                C,
                D,
                U,
                L,
                Preg;
                maxit = maxit,
                tol = tol,
                damping = damping,
                verbose = verbose
            )


        info = (
            method = :litx_regularized_anchor,
            rankA = rA,
            qr_error = qrA.error,
            rankB = infoB.rank,

            regularized = true,
            lambda = λused,

            rankM0 = infoM0.rank,
            M0_fullrank = true,
            sigma_min_M0 = infoM0.σmin,

            iterations = sol.iterations,
            error = sol.error,

            residual_P = sol.residual_P,
            residual_Q = sol.residual_Q
        )


        return sol.P,sol.Q,info
    end


    # ========================================================
    # CASE 2B:
    # M0 still rank deficient
    #
    # User-requested fallback -> ordinary LIT
    # ========================================================

    if verbose

        println(
            "M0 is still rank deficient."
        )

        println(
            "Fallback -> ordinary LIT."
        )
    end


    P,Q =
        LinearTimeIiteration(
            A,
            B,
            C,
            D;
            maxit = maxit,
            tol = tol
        )


    info = (
        method = :lit_fallback,
        rankA = rA,
        qr_error = qrA.error,
        rankB = infoB.rank,

        regularized = true,
        lambda = λused,

        rankM0 = infoM0.rank,
        M0_fullrank = false,
        sigma_min_M0 = infoM0.σmin
    )


    return P,Q,info
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
