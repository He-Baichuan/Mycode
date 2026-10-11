function test_HJB_v4()
%TEST_HJB_V4 Analytic obstacle, QVI, robust saddle, interpolation and LCP signs.
% Requires CompEcon basis functions and lcpsolve. MATLAB-verified 2026-10-10.
run_setup();
fd=struct('maxit',120,'tol',1e-9,'res_tol',1e-7,...
    'Delta',1000,'showiters',false,'lcp_tol',1e-8);
co=struct('maxit',70,'tol',1e-9,'res_tol',1e-7,...
    'showiters',false,'lcp_tol',1e-8);
N=61; F=fundefn('lin',N,0,1); a=funnode(F);
FC=fundefn('cheb',31,0,1); ac=funnode(FC);
% Exact linear complementarity: equation and backend sign conventions.
B=sparse([3,-1;-1,3]); rhs=[0;4]; l=[2;0];
[v,li]=hjb_lcp_solve(B,rhs,l,zeros(2,1),struct('backend','compecon','tol',1e-8));
assert(li.converged && norm(min(v-l,B*v-rhs),inf)<1e-8);
if exist('LCP','file')==2
    [vm,lm]=hjb_lcp_solve(B,rhs,l,zeros(2,1),struct('backend','moll','tol',1e-8));
    assert(lm.converged && norm(v-vm,inf)<1e-6);
end
% Analytic optimal stopping kink.
m=example_model_optstop_v4();
[V,~,rr,id]=HJBsolve(m,F,a,[],fd);
expected=max(10,20*a);
assert(id.converged && norm(V-expected,inf)<1e-7 && norm(rr,inf)<1e-7);
assert(all(id.stop_mask(a<.49)) && ~any(id.stop_mask(a>.51)));
[~,~,W,~,rc,ic]=scsolve_jump(m,FC,ac,[],[],co);
assert(ic.converged && norm(W-max(10,20*ac),inf)<1e-4);
assert(norm(rc,inf)<1e-7);
% Nonlocal regime-switching QVI with known values (17,20).
m=example_model_qvi_v4(); exact=repmat([17,20],N,1);
[V,~,rr,id]=HJBsolve(m,F,a,[],fd);
assert(id.converged && norm(V-exact,inf)<1e-7 && norm(rr(:),inf)<1e-7);
assert(all(id.intervention_mask(:,1)) && ~any(id.intervention_mask(:,2)));
assert(all(id.destination(:,1)==2));
[~,~,W,~,rc,ic]=scsolve_jump(m,FC,ac,[],[],co);
assert(ic.converged && norm(W-repmat([17,20],numel(ac),1),inf)<1e-7);
assert(norm(rc(:),inf)<1e-7);
% Robust stop: 15, 50/3 (rate h=.2 in 2->1 is worst).
m=example_model_robust_stop_v4(); exact=repmat([15,50/3],N,1);
[V,~,rr,id]=HJBsolve(m,F,a,[],fd);
assert(id.converged && norm(V-exact,inf)<1e-6 && norm(rr(:),inf)<1e-7);
assert(all(id.stop_mask(:,1)) && all(abs(id.adversary{2}-.2)<1e-12));
[~,~,W,~,rc,ic]=scsolve_jump(m,FC,ac,[],[],co);
assert(ic.converged && norm(W-repmat([15,50/3],numel(ac),1),inf)<1e-6);
assert(norm(rc(:),inf)<1e-7);
% Generic arbitrary-row interpolation conserves constant functions.
G=hjb_grid(F,a,'fd'); T=[.123;.567;.999];
Y=hjb_value_eval(G,20*ones(N,2),T,[1;2;1],'error');
assert(norm(Y-20,inf)<1e-12);
fprintf('v4 checks PASSED: LCP, optimal stopping, QVI, robust stop, interpolation.\n');
end
