function r=hjb_eval_rho(model,S,j)
N=size(S,1);
if ~isfield(model,'rho'), error('hjb:rho','model.rho required'); end
if isa(model.rho,'function_handle')
    params={}; if isfield(model,'params'), params=model.params; end
    r=model.rho(S,j,params{:});
else
    r=model.rho;
    if isvector(r) && numel(r)>1 && numel(r)~=N, r=r(j); end
    if ismatrix(r) && size(r,2)>1, r=r(:,j); end
end
if isscalar(r), r=repmat(r,N,1); end
if numel(r)~=N || any(~isfinite(r(:))) || any(r(:)<=0)
    error('hjb:rho','Discount must be positive scalar or N-vector per regime');
end
r=r(:);
end
