function [V,x,resid,info,A]=HJBsolve(model,fspace,snodes,V0,options)
%HJBSOLVE Stationary HJB / robust HJB / VI / impulse QVI by finite differences.
% Preserves Moll's SINGLE implicit pseudo-time step per policy update.
% For obstacle-free HJB: [(rho+1/Delta)I-A]Vnew=f+Vold/Delta.
% For stopping/QVI: solve the equivalent frozen-obstacle LCP.
% Economic policy and nature's min choices use shared hjb_stage.
if nargin<3, snodes=[]; end
if nargin<4, V0=[]; end
if nargin<5, options=[]; end
opts=hjb_options(options,'fd');
G=hjb_grid(fspace,snodes,'fd'); R=hjb_regimes(model); N=G.N*R;
if isempty(V0), V=zeros(G.N,R);
elseif numel(V0)==N, V=reshape(V0,G.N,R);
else, error('hjb:V0','Initial V must be N-by-regimes'); end
if any(~isfinite(V(:))), error('hjb:V0','Initial V must be finite'); end
converged=false; distance=zeros(opts.maxit,1); residuals=zeros(opts.maxit,1);
lcpinfo=[];
for iter=1:opts.maxit
    old=V;
    stage=hjb_stage(model,G,old,opts);
    A=stage.A;
    B=spdiags(stage.rho(:)+1/opts.Delta,0,N,N)-A;
    b=stage.f(:)+old(:)/opts.Delta;
    [Gbar,~]=hjb_intervention_operator(model,G,old,opts);
    if any(isfinite(Gbar(:)))
        lopts=struct('backend',opts.lcp_backend,'tol',opts.lcp_tol);
        [solution,lcpinfo]=hjb_lcp_solve(B,b,Gbar(:),old(:),lopts);
    else
        solution=B\b;
    end
    if any(~isfinite(solution(:))), error('hjb:numerics','Nonfinite value solution'); end
    updated=reshape(solution,G.N,R);
    V=(1-opts.qvi_damping)*old+opts.qvi_damping*updated;
    distance(iter)=max(abs(V(:)-old(:)));
    latest=hjb_stage(model,G,V,opts);
    [~,residuals(iter)]=hjb_vi_residual(model,G,latest,V,opts);
    if opts.showiters
        fprintf('FD implicit %3d: dV=%.3e, HJB/VI/QVI-res=%.3e\n',...
            iter,distance(iter),residuals(iter));
    end
    if distance(iter)<opts.tol && (~opts.stop_on_residual || residuals(iter)<opts.res_tol)
        converged=true; break;
    end
end
stage=hjb_stage(model,G,V,opts); A=stage.A;
[resid,rmax,mask]=hjb_vi_residual(model,G,stage,V,opts);
[Gbar,ometa]=hjb_intervention_operator(model,G,V,opts);
x=hjb_pack_policy(stage.policy);
info=struct('method','implicit_fd','converged',converged,...
    'iterations',iter,'error_history',distance(1:iter),...
    'residual_history',residuals(1:iter),'residual_max',rmax,...
    'branches',stage.branch,'grid',{G.nodes},'lcp',lcpinfo,...
    'stop_mask',mask & ometa.kind==1,...
    'intervention_mask',mask & ometa.kind==2,...
    'obstacle_mask',mask,'obstacle',Gbar,'target',ometa.target,...
    'destination',ometa.to,'intervention_choice',ometa.choice,...
    'intervention_action',{ometa.action},'adversary',{stage.adversary},...
    'adversary_branch',stage.adversary_branch,...
    'candidates_per_regime',stage.candidates_count);
if ~converged
    warning('hjb:nonconvergence','FD solver did not meet value + HJB/VI/QVI tolerances');
end
end
