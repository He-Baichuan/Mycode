function [v,info]=hjb_lcp_solve(B,b,lower,v0,options,upper)
%HJB_LCP_SOLVE Standard-sign obstacle LCP, using CompEcon or Moll backends.
% Standard convention: l<=v<=u; F=B*v-b;
%   v=l => F>=0; l<v<u => F=0; v=u => F<=0.
% For lower-bound-only obstacle: 0<=v-l perp F>=0.
% CompEcon lcpsolve uses the OPPOSITE sign, hence (-B,+b).
if nargin<5 || isempty(options), options=struct(); end
if nargin<6 || isempty(upper), upper=inf(size(lower)); end
n=numel(b); b=b(:);
if isscalar(lower), lower=repmat(lower,n,1); else, lower=lower(:); end
if isscalar(upper), upper=repmat(upper,n,1); else, upper=upper(:); end
if nargin<4 || isempty(v0), v0=zeros(n,1); end
v0=v0(:);
if ~isequal(size(B),[n,n]) || numel(lower)~=n || numel(upper)~=n || numel(v0)~=n
    error('hjb:lcp','Inconsistent LCP dimensions');
end
if any(lower>upper) || any(isnan(lower)) || any(isnan(upper)) || any(~isfinite(b))
    error('hjb:lcp','Invalid bounds or RHS');
end
if ~isfield(options,'backend'), options.backend='compecon'; end
if ~isfield(options,'tol'), options.tol=1e-7; end
if all(isinf(lower)&lower<0) && all(isinf(upper)&upper>0)
    v=B\b;
else
    v0=max(lower,min(upper,v0));
    switch char(options.backend)
      case 'compecon'
        if exist('lcpsolve','file')~=2
            error('hjb:lcp','CompEcon lcpsolve.m not on MATLAB path');
        end
        % CompEcon's default smooth transform can have 0/0 Jacobians at
        % degenerate obstacle nodes (x=lower, F=0). Use sparse minmax.
        % Restore caller's CompEcon setting automatically when finished.
        transform=optget('lcpsolve','type','smooth');
        restore=onCleanup(@() optset('lcpsolve','type',transform)); %#ok<NASGU>
        optset('lcpsolve','type','hjb_lcp_minmax');
        [v,~]=lcpsolve(-B,b,lower,upper,v0);
      case 'moll'
        if exist('LCP','file')~=2
            error('hjb:lcp','Yuval Tassa / Moll LCP.m not on MATLAB path');
        end
        v=LCP(B,-b,lower,upper,v0,false);
      otherwise
        error('hjb:lcp','Unknown backend: %s',options.backend);
    end
end
if any(~isfinite(v)), error('hjb:lcp','LCP returned nonfinite values'); end
F=B*v-b;
% Projected normal-map residual: zero iff the two-sided MCP holds.
projected=min(upper,max(lower,v-F));
residual=max(abs(v-projected));
info=struct('backend',options.backend,'residual',residual,...
    'converged',residual<=options.tol,...
    'active_lower',isfinite(lower)&abs(v-lower)<=options.tol,...
    'active_upper',isfinite(upper)&abs(v-upper)<=options.tol,...
    'F',F);
if ~info.converged
    warning('hjb:lcp','LCP backend did not meet complementarity tolerance (%.3g)',residual);
end
end
