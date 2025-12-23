module Numerical_Optim
export GoldenSearch


    function GoldenSearch(f, xl::Float64, xr::Float64; tol = 1e-8)
        p = (sqrt(5.0) - 1.0)/2.0
        b = p*xl + (1.0 -p )*xr
        c = (1.0 - p)*xl + p*xr
        fb = f(b)
        fc = f(c)
        while abs(xr - xl)>tol*max(1.0, abs(b)+abs(c))
            if fb > fc
            xr = c 
            c = b
            fc = fb
            b = p*xl + (1.0 -p )*xr
            fb = f(b)
            else
                xl = b 
                b = c
                fb = fc  
                c = (1.0 - p)*xl + p*xr
                fc = f(c)
            end
        end
        return b 
    end
end
