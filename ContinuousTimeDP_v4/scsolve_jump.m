function [cv,scoord,v,x,resid,info]=scsolve_jump(model,fspace,snodes,v0,x0,options)
%SCSOLVE_JUMP Stationary collocation HJB / robust HJB / VI / impulse QVI.
% Standard full policy evaluation; no pseudo-time iteration.
% For VI/QVI, transform collocation coefficients into nodal values and
% solve node-obstacle LCP with CompEcon lcpsolve (or Moll LCP).
% No guarantee of monotonicity for global Chebyshev bases.
if nargin<3, snodes=[]; end
if nargin<4, v0=[]; end
if nargin<5, x0=[]; end
if nargin<6, options=[]; end
if ~isempty(x0)
    warning('hjb:deprecated','x0 ignored: candidates callback provides policy candidates');
end
opts=hjb_options(options,'collocation');
G=hjb_grid(fspace,snodes,'collocation');
R=hjb_regimes(model); N=G.N*R;
if isempty(v0), V=zeros(G.N,R);
elseif numel(v0)==N, V=reshape(v0,G.N,R);
else, error('hjb:V0','Initial V must be N-by-regimes'); end
if any(~isfinite(V(:))), error('hjb:V0','Initial V must be finite'); end
PB=kron(speye(R),G.P0);
converged=false; distance=zeros(opts.maxit,1); residuals=zeros(opts.maxit,1);
lcpinfo=[];
for iter=1:opts.maxit
    old=V;
    stage=hjb_stage(model,G,old,opts);
    Bc=spdiags(stage.rho(:),0,N,N)*PB-stage.A;
    b=stage.f(:);
    [Gbar,~]=hjb_intervention_operator(model,G,old,opts);
    if any(isfinite(Gbar(:)))
        check_stop_boundary(model,G);
        % Nodal value W=PB*c; complementarity applies to W, NOT c.
        Bvalues=Bc/PB;
        lopts=struct('backend',opts.lcp_backend,'tol',opts.lcp_tol);
        [newV,lcpinfo]=hjb_lcp_solve(Bvalues,b,Gbar(:),old(:),lopts);
        candidate=reshape(newV,G.N,R);
    else
        [Bc,b]=hjb_collocation_neumann(Bc,b,G,R,model);
        candidate=reshape(PB*(Bc\b),G.N,R);
    end
    if any(~isfinite(candidate(:))), error('hjb:numerics','Nonfinite collocation solution'); end
    V=(1-opts.collocation_damping)*old+opts.collocation_damping*candidate;
    distance(iter)=max(abs(V(:)-old(:)));
    latest=hjb_stage(model,G,V,opts);
    [~,residuals(iter)]=hjb_vi_residual(model,G,latest,V,opts);
    if opts.showiters
        fprintf('Collocation PFI %3d: dV=%.3e, HJB/VI/QVI-res=%.3e\n',...
            iter,distance(iter),residuals(iter));
    end
    if distance(iter)<opts.tol && (~opts.stop_on_residual || residuals(iter)<opts.res_tol)
        converged=true; break;
    end
end
cv=G.P0\V;
final=hjb_stage(model,G,V,opts);
[resid,rmax,mask]=hjb_vi_residual(model,G,final,V,opts);
[Gbar,ometa]=hjb_intervention_operator(model,G,V,opts);
scoord=G.xgrid; v=V; x=hjb_pack_policy(final.policy);
info=struct('method','collocation_policy_iteration','converged',converged,...
    'iterations',iter,'error_history',distance(1:iter),...
    'residual_history',residuals(1:iter),'residual_max',rmax,...
    'branches',final.branch,'grid',{G.nodes},'lcp',lcpinfo,...
    'stop_mask',mask & ometa.kind==1,...
    'intervention_mask',mask & ometa.kind==2,...
    'obstacle_mask',mask,'obstacle',Gbar,'target',ometa.target,...
    'destination',ometa.to,'intervention_choice',ometa.choice,...
    'intervention_action',{ometa.action},'adversary',{final.adversary},...
    'adversary_branch',final.adversary_branch,...
    'candidates_per_regime',final.candidates_count);
if ~converged
    warning('hjb:nonconvergence','Collocation solver did not meet value + HJB/VI/QVI tolerances');
end
end
function check_stop_boundary(model,G)
if ~isfield(model,'boundary') || isempty(model.boundary), return; end
bc=model.boundary;
if ischar(bc), bc=repmat({bc},1,G.d); end
for k=1:numel(bc)
    if strcmpi(bc{k},'reflecting')
        error('hjb:collocation_boundary',...
           'Collocation obstacle + reflecting Neumann BC needs augmented active-set boundary equations');
    end
end
end
