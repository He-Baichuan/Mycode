function [residual,rmax,mask]=hjb_vi_residual(model,G,stage,W,opts)
%HJB_VI_RESIDUAL Recompute CURRENT obstacle: full stationary HJB-VI/QVI residual.
% residual = max( obstacle(W)-W, f+A*W-rho.*W ), with -Inf = no obstacle.
N=G.N; R=size(W,2);
state=W;
if strcmp(G.method,'collocation'), state=G.P0\W; end
continuation=reshape(stage.f(:)+stage.A*state(:)-stage.rho(:).*W(:),N,R);
[obstacle,~]=hjb_intervention_operator(model,G,W,opts);
residual=max(obstacle-W,continuation);
if strcmp(G.method,'collocation') && ~any(isfinite(obstacle(:)))
    PB=kron(speye(R),G.P0);
    M=spdiags(stage.rho(:),0,N*R,N*R)*PB-stage.A;
    [M,b]=hjb_collocation_neumann(M,stage.f(:),G,R,model);
    residual=reshape(b-M*state(:),N,R);
end
rmax=max(abs(residual(:)));
mask=isfinite(obstacle) & W-obstacle<=opts.obstacle_tol;
end
