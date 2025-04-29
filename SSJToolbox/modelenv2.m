% model_utils_proc.m: Procedural conversion of the Julia ModelUtils module
% Provides functions to manage model environment, steady states, indexing,
% reshaping, Jacobians, IRFs, nonlinear solves, and optimal policy.

function m = modelenv(par, varFields, steadystate, T)
% MODELENV  Create a model environment struct
%  m = modelenv(par, varFields, steadystate, T)
% Inputs:
%   par        - struct of parameters
%   varFields  - containers.Map mapping keys '0','-1','1',...,'exog' to cell arrays of var names
%   steadystate- struct with fields initial(nx×1), terminal(nx×1), exog(nexog×1)
%   T          - time periods
% Output:
%   m          - struct with fields par, vars, T, ss, nx, nexog, ss_sym
    m.par    = par;
    m.vars   = varFields;
    m.T      = T;
    m.ss     = steadystate;
    m.nx     = numel(steadystate.initial);
    m.nexog  = numel(steadystate.exog);
    if ~isfield(steadystate,'terminal')
        m.ss.terminal = steadystate.initial;
    end
    m.ss_sym = struct();
end

function val = steadystatevalue(s, m, endogexog, which)
% STEADYSTATEVALUE  Get steady-state of variable s
%  val = steadystatevalue(s, m, endogexog, which)
    if nargin<4, which='terminal'; end
    if strcmp(endogexog,'endog')
        vec   = m.ss.(which);
        names = m.vars('0');
    elseif strcmp(endogexog,'exog')
        vec   = m.ss.exog;
        names = m.vars('exog');
    else
        error('endogexog must be ''endog'' or ''exog''.');
    end
    idx = find(strcmp(names,s));
    assert(~isempty(idx), 'Variable not found');
    val = vec(idx);
end

function idx = varindex(s, m, endogexog)
% VARINDEX  Indices of s in long vector X or E
%  idx = varindex(s, m, endogexog)
    if nargin<3, endogexog='endog'; end
    if strcmp(endogexog,'endog')
        names = m.vars('0'); n = m.nx;
    else
        names = m.vars('exog'); n = m.nexog;
    end
    j = find(strcmp(names,s));
    assert(~isempty(j), 'Variable not found');
    idx = (j-1)*m.T + (1:m.T);
end

function [Xss, Ess] = longsteadystate(m)
% LONGSTEADYSTATE  Build long steady-state vectors
    Xss = repmat(m.ss.terminal(:), m.T, 1);
    Ess = repmat(m.ss.exog(:),     m.T, 1);
end

function blocks = wide(x, n)
% WIDE  Split long vector into n blocks of length T
    T = numel(x)/n;
    blocks = mat2cell(x, T*ones(1,n), 1);
end

function xlong = long(blocks)
% LONG  Concatenate cell blocks into a single vector
    xlong = vertcat(blocks{:});
end

function S = contemp(X, m)
% CONTEMP  Struct of contemporaneous values from X
    names  = m.vars('0');
    blocks = wide(X, m.nx);
    for i = 1:m.nx
        S.(names{i}) = blocks{i}(1);
    end
end

function S = lag(X, m, j)
% LAG  Struct of j-period lags
    if nargin<3, j=1; end
    names  = m.vars(num2str(-j));
    blocks = wide(X, m.nx);
    for i = 1:m.nx
        xi = blocks{i};
        S.(names{i}) = [repmat(m.ss.initial(i), j, 1); xi(1:end-j)];
    end
end

function S = lead(X, m, j)
% LEAD  Struct of j-period leads
    if nargin<3, j=1; end
    names  = m.vars(num2str(j));
    blocks = wide(X, m.nx);
    for i = 1:m.nx
        xi = blocks{i};
        S.(names{i}) = [xi(j+1:end); repmat(m.ss.terminal(i), j, 1)];
    end
end

function S = exogenous(E, m)
% EXOGENOUS  Struct of exogenous series from E
    names  = m.vars('exog');
    blocks = wide(E, m.nexog);
    for i = 1:m.nexog
        S.(names{i}) = blocks{i}(1);
    end
end

function J = get_fX(f, m, X, E)
% GET_FX  Numeric Jacobian of f w.r.t. X
    if nargin<4, [X,E]=longsteadystate(m); end
    h = 1e-6; fx0 = f(m,X,E);
    n = numel(X); J = zeros(numel(fx0), n);
    for i = 1:n
        dX = zeros(n,1); dX(i)=h;
        J(:,i) = (f(m,X+dX,E)-fx0)/h;
    end
end

function J = get_fE(f, m, X, E)
% GET_FE  Numeric Jacobian of f w.r.t. E
    if nargin<4, [X,E]=longsteadystate(m); end
    h = 1e-6; fe0 = f(m,X,E);
    n = numel(E); J = zeros(numel(fe0), n);
    for i = 1:n
        dE = zeros(n,1); dE(i)=h;
        J(:,i) = (f(m,X,E+dE)-fe0)/h;
    end
end

function IRF = linearIRFs(f, m, X, E)
% LINEARIRFS  Compute IRFs via -fX\fE
    if nargin<3, [X,E]=longsteadystate(m); end
    FX = get_fX(f,m,X,E);
    FE = get_fE(f,m,X,E);
    IRF = - (FX \ FE);
end

function X = nonlineartransition(m, f, X0, E, maxit)
% NONLINEARTRANSITION  Solve nonlinear equations by Newton
    if nargin<5, maxit=100; end
    X = X0;
    for it = 1:maxit
        res = f(m,X,E);
        FX  = get_fX(f,m,X,E);
        Xold = X;
        X = X - FX \ res;
        if max(abs(X-Xold))<1e-6, return; end
    end
    error('Did not converge');
end

function plotIRFs(IRF, m, irf_horizon, shock, shock_horizon, series)
% PLOTIRFS  Plot IRFs for each endogenous variable
    if nargin<3, irf_horizon=20; end
    if nargin<4, shock=1; end
    if nargin<5, shock_horizon=0; end
    if nargin<6, series=m.vars('0'); end
    % determine shock index
    if isnumeric(shock), sidx=shock;
    else, sidx=find(strcmp(m.vars('exog'),shock)); end
    vec = IRF(:,(sidx-1)*m.T+shock_horizon+1);
    blocks = wide(vec,m.nx);
    figure;
    for i=1:m.nx
        subplot(ceil(m.nx/2),2,i);
        plot(0:irf_horizon, blocks{i}(1:irf_horizon+1));
        title(series{i});
    end
end

function [Xopt, lam] = optimaltransitionpath(m, Ugrad, Uhess, f, X0, E, maxit)
% OPTIMALTRANSITIONPATH  Compute optimal path via KKT
    if nargin<7, maxit=100; end
    X = X0;
    FX = get_fX(f,m,X,E);
    UX = Ugrad(m,X,E);
    lam = FX' \ UX;
    neq = size(FX,1);
    for it=1:maxit
        FX = get_fX(f,m,X,E);
        UX = Ugrad(m,X,E);
        res1 = UX - FX' * lam;
        res2 = f(m,X,E);
        H = numericHessian(@(x) Ugrad(m,x,E)-FX'*lam, X);
        K = [H, -FX'; FX, zeros(neq)];
        upd = -K \ [res1; res2];
        X = X + upd(1:numel(X));
        lam = lam + upd(numel(X)+1:end);
        if max(abs(upd))<1e-7, break; end
    end
    Xopt = X;
end

function [xopt, lam] = optimalLQpolicyExact(Q, R, A, B)
% OPTIMALLQPOLICYEXACT  Solve 0.5 x'Qx + R'x s.t. A + Bx = 0
    K   = [Q, -B'; B, zeros(size(B,1))];
    sol = -K \ [R; A];
    nX  = size(Q,1);
    xopt= sol(1:nX);
    lam = sol(nX+1:end);
end

function H = numericHessian(fun, x)
% NUMERICHESSIAN  Approximate Hessian of vector-valued fun
    n = numel(x); h=1e-5;
    g0 = fun(x);
    H = zeros(n);
    for i=1:n
        ei    = zeros(n,1); ei(i)=h;
        g1    = fun(x+ei);
        H(:,i)= (g1-g0)/h;
    end
end
