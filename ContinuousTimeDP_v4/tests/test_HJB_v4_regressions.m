function test_HJB_v4_regressions()
% Regressions for tensor grids, scalar metadata, boundaries and controls.
run_setup();
opts=struct('showiters',false,'maxit',100,'tol',1e-9,'res_tol',1e-7);
% Solve an analytic constant model with genuine two-dimensional bases.
m=example_model_constant(1);
m.boundary={'state_constraint','state_constraint'};
F=fundefn('lin',[13 9],[0 0],[1 1]); G=hjb_grid(F,[],'fd');
assert(isscalar(G) && iscell(G.xgrid) && numel(G.xgrid)==2);
[V,~,r,info]=HJBsolve(m,F,[],[],opts);
assert(isscalar(info) && iscell(info.grid) && numel(info.grid)==2);
assert(info.converged && max(abs(V-1/m.rho),[],'all')<1e-7 && max(abs(r(:)))<1e-7);
F=fundefn('cheb',[7 6],[0 0],[1 1]); G=hjb_grid(F,[],'collocation');
assert(isscalar(G) && isequal(size(G.P0),[42 42]));
[~,~,V,~,r,info]=scsolve_jump(m,F,[],[],[],opts);
assert(isscalar(info) && iscell(info.grid) && numel(info.grid)==2);
assert(info.converged && max(abs(V-1/m.rho),[],'all')<1e-7 && max(abs(r(:)))<1e-7);
% Nonconstant reflecting solution: check PDE inside and Neumann at edges.
m=struct('rho',.1,'nregimes',1,'boundary','reflecting','candidates',@flow);
F=fundefn('cheb',9,0,1); G=hjb_grid(F,[],'collocation');
[c,~,V,~,r,info]=scsolve_jump(m,F,[],[],[],opts);
assert(info.converged && norm(r,inf)<1e-7);
edge=G.subs(:,1)==1 | G.subs(:,1)==G.n(1);
stage=hjb_stage(m,G,V,hjb_options(opts,'collocation'));
pde=stage.f+stage.A*c-m.rho*V;
assert(max(abs(pde(~edge)))<1e-7);
assert(max(abs(G.P1{1}(edge,:)*c))<1e-7);
assert(max(abs(pde(edge)))>1e-7);
% Independently optimize consumption at the borrowing constraint.
[m,~,~,~,~]=example_model_humancapital(); p=m.params{1};
ellstar=p.alpha*p.delta/(p.rho+p.delta);
hstar=(ellstar^p.alpha*p.theta/p.delta)^(1/(1-p.alpha));
S=[0,hstar]; grad=[.04,.5];
D=struct('forward',grad,'backward',grad,'central',grad,'grad',grad);
cs=m.candidates(S,D,0,1,p); supplied=-inf;
for k=1:numel(cs)
    if cs(k).feasible && cs(k).g(1)>=-1e-11
        supplied=max(supplied,cs(k).f+sum(cs(k).g.*grad,2));
    end
end
[~,objective]=fminbnd(@(ell)-boundary_score(ell,S,grad,p),0,1-1e-10,...
    optimset('TolX',1e-12));
assert(abs(supplied+objective)<1e-8);
fprintf('v4 regressions PASSED: tensor solvers, scalar info, Neumann residual, constrained control.\n');
end

function out=flow(S,~,~,~)
N=size(S,1);
out=struct('x',zeros(N,1),'f',1+S(:,1),'g',zeros(N,1),'sigma',.2*ones(N,1));
end

function H=boundary_score(ell,S,grad,p)
h=S(2); cash=p.w*h*(1-ell); c=min(grad(1)^(-1/p.gamma),cash);
H=c^(1-p.gamma)/(1-p.gamma)+grad(1)*(cash-c)+...
    grad(2)*(p.theta*(ell*h)^p.alpha-p.delta*h);
end
