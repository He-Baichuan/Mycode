function [R,J]=hjb_lcp_minmax(x,a,b,F,Fx)
%HJB_LCP_MINMAX Sparse equivalent of CompEcon minmax.m.
% Keeps Newton Jacobian sparse instead of constructing dense -eye(n).
% Used inside CompEcon lcpsolve via optset('lcpsolve','type',...).
n=numel(x);
if isscalar(a), a=repmat(a,n,1); end
if isscalar(b), b=repmat(b,n,1); end
da=a(:)-x(:); db=b(:)-x(:); F=F(:);
R=min(max(F,da),db);
if nargout>1
    if nargin<5, error('hjb:lcp','Jacobian requested but absent'); end
    J=-speye(n);
    i=find(F>da & F<db);
    if ~isempty(i), J(i,:)=Fx(i,:); end
end
end
