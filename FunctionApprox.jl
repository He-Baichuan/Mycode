module FunctionApprox

export lip, blip, cspline, splint, chebcoef, chebeval1

# ------------------------------- helpers --------------------------------

"""
    _handle_oob(on_extrap, msg)

on_extrap = :error 或 :missing
"""
@inline function _handle_oob(on_extrap::Symbol, msg::String)
    if on_extrap === :missing
        return missing
    elseif on_extrap === :error
        throw(DomainError(msg))
    else
        throw(ArgumentError("on_extrap must be :error or :missing"))
    end
end

@inline function _as_bounds(d)
    a, b = d isa Tuple ? d : (d[1], d[2])
    a < b || throw(ArgumentError("Require lower<upper, got ($a,$b)."))
    return a, b
end

# =============================== LIP ====================================

"""
    lip(xvec, yvec, x0; on_extrap=:error)

一维线性插值：给定严格递增网格 xvec 与函数值 yvec，返回 f(x0) 的线性插值。
x0 可为标量或向量。
"""
function lip(xvec::AbstractVector{<:Real},
             yvec::AbstractVector{<:Real},
             x0; on_extrap::Symbol = :error)

    n = length(xvec)
    n == length(yvec) || throw(DimensionMismatch("xvec and yvec must have same length"))
    n >= 2 || throw(ArgumentError("Need at least 2 grid points"))

    xmin, xmax = xvec[1], xvec[end]

    # scalar
    if x0 isa Real
        x = x0
        if x < xmin || x > xmax
            return _handle_oob(on_extrap, "Input out of grid in lip: x0=$x not in [$xmin,$xmax]")
        end
        if x == xmin
            return yvec[1]
        elseif x == xmax
            return yvec[end]
        else
            j = searchsortedlast(xvec, x)  # largest j with xvec[j] <= x
            # j in 1:(n-1)
            @inbounds begin
                xj, xjp = xvec[j], xvec[j+1]
                yj, yjp = yvec[j], yvec[j+1]
                return yj + (yjp - yj) / (xjp - xj) * (x - xj)
            end
        end
    end

    # vector
    xs = x0::AbstractVector
    m = length(xs)
    out = Vector{Float64}(undef, m)
    @inbounds for k in 1:m
        x = xs[k]
        if x < xmin || x > xmax
            r = _handle_oob(on_extrap, "Input out of grid in lip: x0=$x not in [$xmin,$xmax]")
            if r === missing
                # if user wants missing, promote output to Any
                return Any[ (xi < xmin || xi > xmax) ? missing : lip(xvec, yvec, xi; on_extrap=:error)
                           for xi in xs ]
            end
        end
        if x == xmin
            out[k] = yvec[1]
        elseif x == xmax
            out[k] = yvec[end]
        else
            j = searchsortedlast(xvec, x)
            xj, xjp = xvec[j], xvec[j+1]
            yj, yjp = yvec[j], yvec[j+1]
            out[k] = yj + (yjp - yj) / (xjp - xj) * (x - xj)
        end
    end
    return out
end

# =============================== BLIP ===================================

"""
    blip(xvec, yvec, zmat, x, y; on_extrap=:error)

二维双线性插值：zmat[i,j]=f(xvec[i], yvec[j])。
"""
function blip(xvec::AbstractVector{<:Real},
              yvec::AbstractVector{<:Real},
              zmat::AbstractMatrix{<:Real},
              x::Real, y::Real; on_extrap::Symbol = :error)

    n, m = length(xvec), length(yvec)
    size(zmat, 1) == n || throw(DimensionMismatch("zmat must have size (length(xvec), length(yvec))"))
    size(zmat, 2) == m || throw(DimensionMismatch("zmat must have size (length(xvec), length(yvec))"))

    xmin, xmax = xvec[1], xvec[end]
    ymin, ymax = yvec[1], yvec[end]

    if x < xmin || x > xmax
        return _handle_oob(on_extrap, "x outside of grid in blip: x=$x not in [$xmin,$xmax]")
    end
    if y < ymin || y > ymax
        return _handle_oob(on_extrap, "y outside of grid in blip: y=$y not in [$ymin,$ymax]")
    end

    i = searchsortedlast(xvec, x)
    j = searchsortedlast(yvec, y)

    @inbounds begin
        if i == n && j == m
            return zmat[n, m]
        elseif i == n && j < m
            u = (y - yvec[j]) / (yvec[j+1] - yvec[j])
            return (1 - u) * zmat[n, j] + u * zmat[n, j+1]
        elseif i < n && j == m
            t = (x - xvec[i]) / (xvec[i+1] - xvec[i])
            return t * zmat[i+1, m] + (1 - t) * zmat[i, m]
        else
            t = (x - xvec[i]) / (xvec[i+1] - xvec[i])
            u = (y - yvec[j]) / (yvec[j+1] - yvec[j])
            return (1 - t) * (1 - u) * zmat[i, j] +
                   t       * (1 - u) * zmat[i+1, j] +
                   t       * u       * zmat[i+1, j+1] +
                   (1 - t) * u       * zmat[i, j+1]
        end
    end
end

# =============================== CSpline ================================

"""
    cspline(x, y, cmethod, yp)

计算三次样条所需的二阶导 y2（对应 Numerical Recipes 的 spline）。
- cmethod=1: natural spline, y2[1]=y2[n]=0
- cmethod=2: 端点一阶导用割线（secant）近似
- cmethod=3: 端点一阶导由 yp=(yp1, ypn) 给定（或长度2向量）
"""
function cspline(x::AbstractVector{<:Real},
                 y::AbstractVector{<:Real},
                 cmethod::Integer;
                 yp = (0.0, 0.0))

    n = length(x)
    n == length(y) || throw(DimensionMismatch("x and y must have same length"))
    n >= 2 || throw(ArgumentError("Need at least 2 points"))
    # assume x strictly increasing as in GAUSS code

    y2 = zeros(Float64, n)
    u  = zeros(Float64, n)

    qn = 0.0
    un = 0.0

    if cmethod == 1
        y2[1] = 0.0
        u[1]  = 0.0
        qn    = 0.0
        un    = 0.0
    elseif cmethod == 2
        yp1 = (y[2] - y[1]) / (x[2] - x[1])
        ypn = (y[n] - y[n-1]) / (x[n] - x[n-1])
        y2[1] = -0.5
        u[1]  = (3.0 / (x[2] - x[1])) * ((y[2] - y[1]) / (x[2] - x[1]) - yp1)
        qn    = 0.5
        un    = (3.0 / (x[n] - x[n-1])) * (ypn - (y[n] - y[n-1]) / (x[n] - x[n-1]))
    elseif cmethod == 3
        yp1, ypn = yp isa Tuple ? yp : (yp[1], yp[2])
        y2[1] = -0.5
        u[1]  = (3.0 / (x[2] - x[1])) * ((y[2] - y[1]) / (x[2] - x[1]) - yp1)
        qn    = 0.5
        un    = (3.0 / (x[n] - x[n-1])) * (ypn - (y[n] - y[n-1]) / (x[n] - x[n-1]))
    else
        throw(ArgumentError("cmethod must be 1,2,3"))
    end

    @inbounds for i in 2:(n-1)
        sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1])
        p   = sig * y2[i-1] + 2.0
        y2[i] = (sig - 1.0) / p
        u[i]  = (6.0 * ((y[i+1] - y[i]) / (x[i+1] - x[i]) -
                        (y[i]   - y[i-1]) / (x[i]   - x[i-1])) / (x[i+1] - x[i-1])
                 - sig * u[i-1]) / p
    end

    @inbounds begin
        y2[n] = (un - qn * u[n-1]) / (qn * y2[n-1] + 1.0)
        for k in (n-1):-1:1
            y2[k] = y2[k] * y2[k+1] + u[k]
        end
    end

    return y2
end

# =============================== Splint ================================

"""
    splint(x, y, y2, gmethod, x0; on_extrap=:error)

用三次样条在点 x0 处插值（需先由 cspline 得到 y2）。
- gmethod=1: 用 searchsortedlast 定位 bracket（等价于 GAUSS 的 sumc(x.<=x0)）
- gmethod=2: 手写二分（基本等价）
x0 可为标量或向量。
"""
function splint(x::AbstractVector{<:Real},
                y::AbstractVector{<:Real},
                y2::AbstractVector{<:Real},
                gmethod::Integer,
                x0; on_extrap::Symbol = :error)

    n = length(x)
    n == length(y) == length(y2) || throw(DimensionMismatch("x,y,y2 must have same length"))
    n >= 2 || throw(ArgumentError("Need at least 2 points"))
    xmin, xmax = x[1], x[end]

    # locate bracket
    @inline function _bracket(xi::Real)
        if gmethod == 1
            klo = searchsortedlast(x, xi)
            klo = min(max(klo, 1), n-1)
            return klo, klo+1
        elseif gmethod == 2
            klo = 1
            khi = n
            while (khi - klo) > 1
                k = (khi + klo) >>> 1
                if x[k] > xi
                    khi = k
                else
                    klo = k
                end
            end
            return klo, khi
        else
            throw(ArgumentError("gmethod must be 1 or 2"))
        end
    end

    # scalar
    if x0 isa Real
        xi = x0
        if xi < xmin || xi > xmax
            return _handle_oob(on_extrap, "x0 out of grid in splint: x0=$xi not in [$xmin,$xmax]")
        end
        if xi == xmin
            return y[1]
        elseif xi == xmax
            return y[end]
        end
        klo, khi = _bracket(xi)
        @inbounds begin
            h = x[khi] - x[klo]
            a = (x[khi] - xi) / h
            b = (xi - x[klo]) / h
            return a*y[klo] + b*y[khi] +
                   ((a^3 - a) * y2[klo] + (b^3 - b) * y2[khi]) * (h^2) / 6.0
        end
    end

    # vector
    xs = x0::AbstractVector
    m = length(xs)
    out = Vector{Float64}(undef, m)
    @inbounds for i in 1:m
        xi = xs[i]
        if xi < xmin || xi > xmax
            r = _handle_oob(on_extrap, "x0 out of grid in splint: x0=$xi not in [$xmin,$xmax]")
            if r === missing
                return Any[ (xj < xmin || xj > xmax) ? missing : splint(x,y,y2,gmethod,xj; on_extrap=:error)
                           for xj in xs ]
            end
        end
        if xi == xmin
            out[i] = y[1]
        elseif xi == xmax
            out[i] = y[end]
        else
            klo, khi = _bracket(xi)
            h = x[khi] - x[klo]
            a = (x[khi] - xi) / h
            b = (xi - x[klo]) / h
            out[i] = a*y[klo] + b*y[khi] +
                     ((a^3 - a) * y2[klo] + (b^3 - b) * y2[khi]) * (h^2) / 6.0
        end
    end
    return out
end

# ============================== ChebCoef ===============================

"""
    chebcoef(f, n, m, d)

Chebyshev 回归（与 GAUSS 代码同节点与同系数公式）：
- n: 系数个数（对应 T_0,...,T_{n-1}）
- m: 采样点数（若 m<n 自动设为 n）
- d: (a,b) 或长度2向量
返回 alpha 向量（长度 n）。
"""
function chebcoef(f::Function, n::Integer, m::Integer, d)
    n >= 1 || throw(ArgumentError("n must be >=1"))
    m < n && (m = n)
    a, b = _as_bounds(d)

    k = collect(1:m)
    xbar = cos.(((2 .* k .- 1) ./ (2m)) .* pi)               # zeros on [-1,1]
    zbar = (xbar .+ 1) .* (b - a) .* 0.5 .+ a                # map to [a,b]
    ybar = f.(zbar)

    alpha = zeros(Float64, n)
    alpha[1] = (1 / m) * sum(ybar)
    if n >= 2
        t1 = xbar
        alpha[2] = (2 / m) * dot(ybar, t1)
    end
    if n >= 3
        t0 = ones(Float64, m)
        t1 = xbar
        t2 = 2 .* xbar .* t1 .- t0
        for i in 3:n
            alpha[i] = (2 / m) * dot(ybar, t2)
            t0, t1 = t1, t2
            t2 = 2 .* xbar .* t1 .- t0
        end
    end
    return alpha
end

# ============================== ChebEval1 ==============================

"""
    chebeval1(alpha, z, d; method=:clenshaw)

在区间 d=(a,b) 上，对点 z 评价 Chebyshev 多项式
    Σ_{ℓ=0}^{n-1} alpha[ℓ+1] T_ℓ(x),  x = 2(z-a)/(b-a)-1.
- method=:clenshaw（推荐，稳定且省内存）
- method=:matrix（对应 GAUSS 的显式构造 T 矩阵）
z 可为标量或向量。
"""
function chebeval1(alpha::AbstractVector{<:Real}, z, d; method::Symbol = :clenshaw)
    a, b = _as_bounds(d)
    n = length(alpha)

    # map
    mapx(zz) = (2 * (zz - a) / (b - a)) - 1

    # scalar
    if z isa Real
        x = mapx(z)
        if method === :matrix
            # build T0..T_{n-1}
            if n == 1
                return float(alpha[1])
            elseif n == 2
                return float(alpha[1] + alpha[2]*x)
            else
                t0 = 1.0
                t1 = x
                acc = alpha[1]*t0 + alpha[2]*t1
                for k in 2:(n-1)
                    t2 = 2x*t1 - t0
                    acc += alpha[k+1]*t2
                    t0, t1 = t1, t2
                end
                return float(acc)
            end
        elseif method === :clenshaw
            if n == 1
                return float(alpha[1])
            end
            bkp1 = 0.0
            bkp2 = 0.0
            @inbounds for j in (n-1):-1:1
                bj  = 2x*bkp1 - bkp2 + alpha[j+1]
                bkp2 = bkp1
                bkp1 = bj
            end
            return float(x*bkp1 - bkp2 + alpha[1])
        else
            throw(ArgumentError("method must be :clenshaw or :matrix"))
        end
    end

    # vector
    zs = z::AbstractVector
    m = length(zs)
    out = Vector{Float64}(undef, m)
    @inbounds for i in 1:m
        out[i] = chebeval1(alpha, zs[i], (a,b); method=method)
    end
    return out
end

end # module
