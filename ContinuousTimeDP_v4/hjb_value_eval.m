function y=hjb_value_eval(G,W,T,to,outside)
%HJB_VALUE_EVAL Evaluate V_j at arbitrary target points T for FD/collocation.
% to: scalar or number-of-targets vector of destination regime IDs.
if nargin<5 || isempty(outside), outside='error'; end
if size(T,2)~=G.d || any(~isfinite(T(:)))
    error('hjb:target','Target must be finite m-by-d');
end
m=size(T,1); R=size(W,2);
if isscalar(to), to=repmat(to,m,1); else, to=to(:); end
if numel(to)~=m || any(to<1 | to>R | to~=round(to))
    error('hjb:target','Invalid destination regime');
end
if strcmp(G.method,'fd')
    P=hjb_fd_interp(G,T,outside);
    y=zeros(m,1);
    for j=1:R
        ix=find(to==j);
        if ~isempty(ix), y(ix)=P(ix,:)*W(:,j); end
    end
else
    for k=1:G.d
        left=G.nodes{k}(1); right=G.nodes{k}(end);
        if strcmpi(outside,'clip')
            T(:,k)=min(max(T(:,k),left),right);
        elseif any(T(:,k)<left-1e-12 | T(:,k)>right+1e-12)
            error('hjb:target','Target outside collocation domain');
        end
    end
    P=funbas(G.fspace,T,zeros(1,G.d));
    C=G.P0\W; y=zeros(m,1);
    for j=1:R
        ix=find(to==j);
        if ~isempty(ix), y(ix)=P(ix,:)*C(:,j); end
    end
end
end
