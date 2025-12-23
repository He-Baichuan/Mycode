module FixvMN1
export fixvmn1
using LinearAlgebra

"""
crit = (status, resid_inf, relchange, obj, iter)

status:
  0 : 正常收敛
  1 : 函数评估/雅可比/线性方程失败（NaN/Inf/异常/奇异）
  2 : 线搜索无法进一步降低目标（NRStep 返回 -1）
  3 : 超过最大迭代次数
"""
Base.@kwdef struct MNRCrit
    status::Int = 0
    resid_inf::Float64 = Inf
    relchange::Float64 = Inf
    obj::Float64 = Inf
    iter::Int = 0
end

# ---------- utilities ----------
_isbad(v) = any(x -> !isfinite(x), v)
_isbad(A::AbstractMatrix) = any(x -> !isfinite(x), A)

# safe evaluation: return nothing if throw/NaN/Inf
function _safeeval(f, x::AbstractVector)
    fx = try
        f(x)
    catch
        return nothing
    end
    return _isbad(fx) ? nothing : fx
end

# Jacobian health (fast for small n): returns (σmin, cond_est)
function _jac_health(J::AbstractMatrix)
    sv = svdvals(Matrix(J))
    σmin = sv[end]
    σmax = sv[1]
    cond_est = σmax / max(σmin, eps(Float64))
    return σmin, cond_est
end
# ParTest(x, dx, typx): maximum relative change
function par_test(x::AbstractVector, dx::AbstractVector, typx::AbstractVector)
    @assert length(x) == length(dx) == length(typx)
    crit = 0.0
    @inbounds for i in eachindex(x)
        denom = max(abs(x[i]), abs(typx[i]))
        crit = max(crit, abs(dx[i]) / denom)
    end
    return crit
end

# CDJac: central difference Jacobian
function cdjac(f, x::AbstractVector{<:Real}; fx0=nothing)
    x0 = Vector{Float64}(x)
    f0 = fx0 === nothing ? f(x0) : fx0
    n = length(f0)
    m = length(x0)
    J = zeros(Float64, n, m)

    epsh = eps(Float64)^(1/3)
    x1 = copy(x0)
    x2 = copy(x0)

    @inbounds for i in 1:m
        # h = sign(x_i) * eps^(1/3) * max(|x_i|, 1)
        h = (x0[i] < 0 ? -1.0 : 1.0) * epsh * max(abs(x0[i]), 1.0)

        temp = x0[i]
        x1[i] = temp + h
        h = x1[i] - temp          # Dennis–Schnabel trick
        x2[i] = temp - h

        f1 = _safeeval(f, x1)
        f2 = _safeeval(f, x2)
        if f1 === nothing || f2 === nothing
            fill!(J, NaN)
            return J
        end


        @inbounds for j in 1:n
            J[j, i] = (f1[j] - f2[j]) / (2.0 * h)
        end

        x1[i] = temp
        x2[i] = temp
    end

    return J
end

# LSolve via QR (as in Toolbox.src)
lsolve_qr(A, b) = (qr(A) \ b)

# Objective φ(x) = 1/2 ||f(x)||^2
phi_from_fx(fx::AbstractVector) = 0.5 * dot(fx, fx)

"""
NRStep: line search on φ(x)=1/2||f(x)||^2 with Armijo-like condition.
Returns step s ∈ (0,1] or -1.0 if no further progress possible.
"""
function nr_step(x0::AbstractVector, dx0::AbstractVector, dg::AbstractVector, f;
                 smult=1.0e-4, smin=0.1, smax=0.5, stol=1.0e-11)

    s1 = 1.0
    fx0 = _safeeval(f, x0)
    fx0 === nothing && return -1.0
    #f(x0)
    g0 = phi_from_fx(fx0)

    dgdx = dot(dg, dx0)  # ∇φ(x0)' dx0
    fx1 = _safeeval(f, x0 .+ dx0)
    g1 = fx1 === nothing ? Inf : phi_from_fx(fx1)
    #fx1 = f(x0 .+ dx0)
    #g1 = phi_from_fx(fx1)

    # accept full step if sufficient decrease
    if g1 <= g0 + smult * dgdx
        return s1
    end

    # quadratic step
    denom = 2.0 * (g1 - g0 - dgdx)
    s = -dgdx / denom
    s = clamp(s, smin, smax)
    fx2 = _safeeval(f, x0 .+ s .* dx0)
    g2 = fx2 === nothing ? Inf : phi_from_fx(fx2)
    #fx2 = f(x0 .+ s .* dx0)
    #g2 = phi_from_fx(fx2)
    s2 = s

    # cubic backtracking loop
    while g2 > (g0 + smult * s2 * dgdx)
        a11 = 1.0 / (s2^2);   a12 = -1.0 / (s1^2)
        a21 = -s1 / (s2^2);   a22 =  s2 / (s1^2)

        b1 = g2 - s2 * dgdx - g0
        b2 = g1 - s1 * dgdx - g0

        # ab = (amat*bvec)/(s2-s1)
        den = (s2 - s1)
        a = (a11*b1 + a12*b2) / den
        b = (a21*b1 + a22*b2) / den

        if a == 0.0
            s = -dgdx / (2.0 * b)
        else
            disc = b^2 - 3.0 * a * dgdx
            if disc < 0.0
                s = s2 * smax
            elseif b <= 0.0
                s = (-b + sqrt(disc)) / (3.0 * a)
            else
                s = -dgdx / (b + sqrt(disc))
            end
        end

        s = min(s, s2 * smax)

        tol = norm(s .* dx0) / (1.0 + norm(x0))
        if tol < stol
            return -1.0
        end

        # shift
        s1 = s2
        s2 = s
        g1 = g2
        fx2 = _safeeval(f, x0 .+ s2 .* dx0)
        g2 = fx2 === nothing ? Inf : phi_from_fx(fx2)
        #fx2 = f(x0 .+ s2 .* dx0)
        #g2 = phi_from_fx(fx2)
    end

    return s2
end

"""
FixvMN1 solver.

Arguments:
  x0::Vector          initial guess
  f                  function f(x)::Vector

Keywords:
  maxit=5000, stopc=1e-8
  global_search=true     (use NRStep)
  use_qr=false           (QR solve for Newton step)
  verbose=false

Returns:
  x, crit::MNRCrit
"""
function fixvmn1(x0::AbstractVector, f;
                 maxit::Int=5000,
                 stopc::Float64=1e-8,
                 global_search::Bool=true,
                 use_qr::Bool=false,
                 verbose::Bool=false)

    x = Vector{Float64}(x0)
    crit = MNRCrit(status=0, resid_inf=Inf, relchange=Inf, obj=Inf, iter=0)

    critold_obj = 2.0

    while !(crit.resid_inf < stopc || crit.iter > maxit)

        # Evaluate f(x)
        fx = try
            f(x)
        catch
            crit = MNRCrit(status=1, resid_inf=crit.resid_inf, relchange=crit.relchange,
                           obj=crit.obj, iter=crit.iter)
            return x, crit
        end
        if _isbad(fx)
            crit = MNRCrit(status=1, resid_inf=crit.resid_inf, relchange=crit.relchange,
                           obj=crit.obj, iter=crit.iter)
            return x, crit
        end
        resid0 = maximum(abs.(fx))
        obj0   = phi_from_fx(fx)
        if verbose
            println("iter=$(crit.iter)  ||f(x)||∞=", resid0, "  φ(x)=", obj0)
        end

        # optional: reject obviously insane starting points
        if resid0 > 1e12
            crit = MNRCrit(status=1, resid_inf=resid0, relchange=Inf, obj=obj0, iter=crit.iter)
            return x, crit
        end
        # Jacobian
        J = cdjac(f, x; fx0=fx)
        if _isbad(J)
            crit = MNRCrit(status=1, resid_inf=crit.resid_inf, relchange=crit.relchange,
                           obj=crit.obj, iter=crit.iter)
            return x, crit
        end

        # gradient of φ(x)=1/2||f||^2 : ∇φ = J' f
        dg = global_search ? (transpose(J) * fx) : zeros(Float64, length(x))
        #=
        # Newton step: J * dx = -f
        dx = try
            use_qr ? lsolve_qr(J, -fx) : (J \ (-fx))
        catch
            crit = MNRCrit(status=1, resid_inf=crit.resid_inf, relchange=crit.relchange,
                           obj=crit.obj, iter=crit.iter)
            return x, crit
        end
        =#
        # Newton step: J * dx = -f  (with singular/ill-conditioned fallback)
        dx = try
            σmin, condJ = _jac_health(J)

            if verbose
                println("    σmin(J)=", σmin, "  cond_est(J)=", condJ)
            end

            if σmin < 1e-12 || condJ > 1e12
                # fallback 1: pivoted QR (robust for rank-deficient / ill-conditioned)
                -(qr(J, Val(true)) \ fx)
                # (optional) fallback 2 (LM), comment out QR above and use:
                # λ = 1e-8 * (svdvals(J)[1]^2)
                # -((J'J + λ*I) \ (J' * fx))
            else
                use_qr ? lsolve_qr(J, -fx) : (J \ (-fx))
            end
        catch
            crit = MNRCrit(status=1, resid_inf=crit.resid_inf, relchange=crit.relchange,
                        obj=crit.obj, iter=crit.iter)
            return x, crit
        end
        # Step1: avoid non-finite evaluations
        step1 = 1.0
        while true
            fx_try = try
                f(x .+ step1 .* dx)
            catch
                nothing
            end
            if fx_try !== nothing && !_isbad(fx_try)
                break
            end
            step1 *= 0.75
            if step1 < 1.0e-16
                crit = MNRCrit(status=1, resid_inf=crit.resid_inf, relchange=crit.relchange,
                               obj=crit.obj, iter=crit.iter)
                return x, crit
            end
        end
        dx .= step1 .* dx

        # Step2: global line search (optional)
        step2 = global_search ? nr_step(x, dx, dg, f) : 1.0
        if step2 < 0.0
            crit = MNRCrit(status=2, resid_inf=crit.resid_inf, relchange=crit.relchange,
                           obj=crit.obj, iter=crit.iter)
            return x, crit
        end

        x_new = x .+ step2 .* dx
        fx_new = f(x_new)

        resid_inf = maximum(abs.(fx_new))
        relchange = par_test(x, step2 .* dx, ones(length(x)))
        obj = phi_from_fx(fx_new)

        if verbose
            println("iter=$(crit.iter)  ||f||_∞=$resid_inf  relchange=$relchange  obj=$obj  step1=$step1  step2=$step2")
            if global_search
                println("    decrease? ", obj < critold_obj)
            end
        end

        critold_obj = obj
        x = x_new
        crit = MNRCrit(status=0, resid_inf=resid_inf, relchange=relchange, obj=obj, iter=crit.iter + 1)
    end

    if crit.iter >= maxit
        crit = MNRCrit(status=3, resid_inf=crit.resid_inf, relchange=crit.relchange, obj=crit.obj, iter=crit.iter)
    end
    return x, crit
end

end # module
