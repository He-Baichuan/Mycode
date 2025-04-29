function m = ModelEnv(par, vars, steadystate, T)
            % Constructor: initialize ModelUtils environment
            % varFields: struct with fields '0','-1','1',...,'exog' mapping to cell arrays of names
            m = ModelUtils;
            m.par = par;
        
            m.vars.endo = vars.endo;
            m.vars.exog = vars.exog;
           
            m.T = T;
            m.nx = numel(steadystate.initial);
            m.nexog = numel(steadystate.exog);

            if ~isfield(steadystate,'terminal')
                steadystate.terminal = steadystate.initial;
            end
            m.ss = steadystate;

            m.idx_now  = containers.Map(vars.endo,  1:m.nx);
            m.idx_exog = containers.Map(vars.exog,1:m.nexog);
        end