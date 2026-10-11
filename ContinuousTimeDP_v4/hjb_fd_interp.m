function P=hjb_fd_interp(G,T,outside)
%HJB_FD_INTERP Multilinear nonnegative interpolation, m queries by N nodes.
% Handles nonuniform axes and exactly-on-boundary targets.
if nargin<3 || isempty(outside), outside='error'; end
N=G.N; d=G.d; m=size(T,1);
if size(T,2)~=d, error('hjb:jump','Jump targets must be N-by-d'); end
if any(~isfinite(T(:))), error('hjb:jump','Jump target contains NaN/Inf'); end
lo=zeros(m,d); frac=zeros(m,d);
for k=1:d
    nodes=G.nodes{k}; x=T(:,k);
    if strcmpi(outside,'clip')
        x=max(nodes(1),min(nodes(end),x));
    elseif any(x<nodes(1)-1e-12 | x>nodes(end)+1e-12)
        error('hjb:jump','Jump target outside grid; specify jump_outside=clip if intentional');
    end
    x=max(nodes(1),min(nodes(end),x));
    raw=interp1(nodes,(1:G.n(k))',x,'linear');
    idx=min(G.n(k)-1,max(1,floor(raw)));
    lo(:,k)=idx;
    frac(:,k)=(x-nodes(idx))./(nodes(idx+1)-nodes(idx));
    frac(:,k)=max(0,min(1,frac(:,k)));
end
rows=cell(2^d,1); cols=rows; vals=rows;
for bits=0:2^d-1
    dest=ones(m,1); wt=ones(m,1);
    for k=1:d
        high=bitget(bits,k);
        dest=dest+(lo(:,k)+high-1)*G.stride(k);
        if high, wt=wt.*frac(:,k); else, wt=wt.*(1-frac(:,k)); end
    end
    rows{bits+1}=(1:m)'; cols{bits+1}=dest; vals{bits+1}=wt;
end
P=sparse(vertcat(rows{:}),vertcat(cols{:}),vertcat(vals{:}),m,N);
end
