function test_HJBkit()
%TEST_HJBKIT Independent analytic identities and operator invariants.
% Requires MATLAB+CompEcon. No claims without running this script.
run_setup();
fd=struct('showiters',false,'maxit',100,'Delta',1000,...
  'tol',1e-9,'res_tol',2e-8);
co=struct('showiters',false,'maxit',60,'tol',1e-9,...
  'res_tol',2e-8);
passed=0;

% 1. Two-state switching, constant flows and nonzero inward drift.
mod1=example_model_constant(2);
F=fundefn('lin',61,0,1); a=funnode(F);
[V,x,r,info,A]=HJBsolve(mod1,F,a,zeros(61,2),fd); %#ok<ASGLU>
exact=(mod1.rho*eye(2)-mod1.Q)\[1;2];
assert(info.converged && max(abs(V-repmat(exact',61,1)),[],'all')<1e-7);
gs=hjb_check_generator(A,1e-9,false); assert(gs.valid);
assert(max(abs(r(:)))<2e-8); passed=passed+1;
FC=fundefn('cheb',18,0,1); ac=funnode(FC);
[C,~,W,~,rc,info]=scsolve_jump(mod1,FC,ac,zeros(18,2),[],co); %#ok<ASGLU>
assert(info.converged && max(abs(W-repmat(exact',18,1)),[],'all')<1e-7);
assert(max(abs(rc(:)))<2e-8); passed=passed+1;

% 2. Nonlocal Poisson reset, analytical affine value.
mod2=example_model_reset(); lam=.4; rho=.1; target=.5;
F=fundefn('lin',55,0,1); a=funnode(F);
[V,~,rr,info,A]=HJBsolve(mod2,F,a,zeros(55,1),fd);
exact=(a+lam*target/rho)/(rho+lam);
assert(info.converged && max(abs(V-exact))<1e-7);
gs=hjb_check_generator(A,1e-9,false); assert(gs.valid); passed=passed+1;
FC=fundefn('cheb',15,0,1); ac=funnode(FC);
[C,~,V,~,rr,info]=scsolve_jump(mod2,FC,ac,zeros(15,1),[],co); %#ok<ASGLU>
exact=(ac+lam*target/rho)/(rho+lam);
assert(info.converged && max(abs(V-exact))<1e-7); passed=passed+1;

% 3. Discrete branch global maximization.
mod3=example_model_branches();
F=fundefn('lin',40,0,1); a=funnode(F);
[V,x,~,info]=HJBsolve(mod3,F,a,zeros(40,1),fd);
assert(info.converged && max(abs(V-20))<1e-7 && all(x(:)==2)); passed=passed+1;
FC=fundefn('cheb',13,0,1); ac=funnode(FC);
[~,~,V,x,~,info]=scsolve_jump(mod3,FC,ac,zeros(13,1),[],co);
assert(info.converged && max(abs(V-20))<1e-7 && all(x(:)==2)); passed=passed+1;

% 4. Reflecting diffusions; constant solution and generator check.
mod4=example_model_diffusion();
F=fundefn('lin',45,0,1); a=funnode(F);
[V,~,rr,info,A]=HJBsolve(mod4,F,a,zeros(45,1),fd);
assert(info.converged && max(abs(V-20))<1e-7);
gs=hjb_check_generator(A,1e-9,false); assert(gs.valid); passed=passed+1;
FC=fundefn('cheb',17,0,1); ac=funnode(FC);
[~,~,V,~,rr,info]=scsolve_jump(mod4,FC,ac,zeros(17,1),[],co);
assert(info.converged && max(abs(V-20))<1e-6); passed=passed+1;

% 5. Two-dimensional NONUNIFORM FD grid.
F=fundefn('lin',[17,11],[0 0],[1 1]);
nd={linspace(0,1,17)'.^1.4,linspace(0,1,11)'.^1.8};
G=hjb_grid(F,nd,'fd');
mu=[0.2*(0.5-G.S(:,1)),0.2*(0.5-G.S(:,2))];
L=hjb_fd_operator(G,mu,[],{'state_constraint','state_constraint'},hjb_options([],'fd'));
gs=hjb_check_generator(L,1e-9,false); assert(gs.valid);
% Standalone implicit system with known constant steady payoff.
rho=.1; val=(rho*speye(G.N)-L)\ones(G.N,1);
assert(max(abs(val-1/rho))<1e-8); passed=passed+1;

% 6. Correlated 2D reflecting diffusions, positive-type covariance.
F=fundefn('lin',[12,14],[0 0],[1 1]); G=hjb_grid(F,[], 'fd');
N=G.N; mu=zeros(N,2); sig=zeros(N,2,2);
sig(:,1,1)=.2; sig(:,2,1)=.03;
sig(:,2,2)=sqrt(.04-.03^2);
L=hjb_fd_operator(G,mu,sig,{'reflecting','reflecting'},hjb_options([],'fd'));
gs=hjb_check_generator(L,1e-9,false); assert(gs.valid);
assert(max(abs((.1*speye(N)-L)\ones(N,1)-10))<1e-8); passed=passed+1;

% 7. Stylized kink / nonconvex branch switch (not Moll's Skiba model).
mod7=example_model_skiba_branches();
F=fundefn('lin',61,0,1); a=funnode(F);
[V,x,rr,info]=HJBsolve(mod7,F,a,zeros(61,1),fd);
truth=max(1+.8*a,1.8-.8*a)/.1;
assert(info.converged && max(abs(V-truth))<1e-7);
assert(all(x(a<.5)==2) && all(x(a>.5)==1)); passed=passed+1;

% 8. Genuine Moll/Kaplan Skiba specification: candidate generation and
% zero-outflow boundaries; this is not a full convergence proof.
[ms,kmin,kmax,~]=example_model_skiba_moll();
F=fundefn('lin',51,kmin,kmax); kk=funnode(F); G=hjb_grid(F,kk,'fd');
p=ms.params{1}; guess=((kk.^p.alpha).^(1-p.s)/(1-p.s))/p.rho;
stage=hjb_stage(ms,G,guess,hjb_options([],'fd'));
gs=hjb_check_generator(stage.A,1e-9,false); assert(gs.valid); passed=passed+1;

% 9. Human-capital joint-control generator test on a modest 2D grid.
[mh,amin,amax,hmin,hmax]=example_model_humancapital();
F=fundefn('lin',[13 9],[amin hmin],[amax hmax]); G=hjb_grid(F,[], 'fd');
p=mh.params{1}; av=G.S(:,1); hv=G.S(:,2);
ell=min(1,max(0,((p.delta/p.theta).*hv.^(1-p.alpha)).^(1/p.alpha)));
c0=p.r*av+p.w*hv.*(1-ell);
v0=c0.^(1-p.gamma)/(1-p.gamma)/mh.rho;
stage=hjb_stage(mh,G,v0,hjb_options([],'fd'));
gs=hjb_check_generator(stage.A,1e-9,false); assert(gs.valid);
assert(size(stage.policy{1},2)==2); passed=passed+1;

fprintf('ALL %d MATLAB verification cases PASSED.\n',passed);
end
