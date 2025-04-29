function IRF = linearIRFs(f, m, X, E)
            if nargin<3
                [X,E] = longsteadystate(m); 
            end
            [FX, FE] = get_jacobian(f, m, X, E);
            IRF =-( sparse(FX) \ FE );
        end