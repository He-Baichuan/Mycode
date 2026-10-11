function G=hjb_grid(fspace,snodes,method)
%HJB_GRID Coordinate ordering: first dimension varies fastest (NDGRID).
if nargin<3, method='fd'; end
if nargin<2 || isempty(snodes), snodes=funnode(fspace); end
if ~iscell(snodes), snodes={snodes}; end
d=numel(snodes); n=zeros(1,d);
for k=1:d
    snodes{k}=snodes{k}(:);
    n(k)=numel(snodes{k});
    if n(k)<3 || any(diff(snodes{k})<=0)
        error('hjb:grid','Each coordinate requires >=3 strictly increasing nodes');
    end
end
if d~=fspace.d || any(n~=reshape(fspace.n,1,[]))
    error('hjb:grid','Grid does not match fspace.n');
end
if strcmpi(method,'fd') && ~all(strcmp(fspace.bastype,'lin'))
    error('hjb:grid','Finite-difference grid must use linear CompEcon family');
end
if d==1
    S=snodes{1}; in=snodes{1};
else
    S=gridmake(snodes); in=snodes;
end
N=size(S,1);
subs=zeros(N,d);
stride=[1 cumprod(n(1:end-1))];
for k=1:d
    subs(:,k)=mod(floor((0:N-1)'/stride(k)),n(k))+1;
end
G=struct('n',n,'d',d,'N',N,'nodes',{snodes},'S',S,...
    'subs',subs,'stride',stride,'fspace',fspace,'xgrid',{in},'method',method);
if strcmpi(method,'collocation')
    G.P0=funbas(fspace,in,zeros(1,d));
    G.P1=cell(1,d); G.P2=cell(d,d);
    for k=1:d
        deriv=zeros(1,d); deriv(k)=1;
        G.P1{k}=funbas(fspace,in,deriv);
        for l=k:d
            deriv2=deriv; deriv2(l)=deriv2(l)+1;
            G.P2{k,l}=funbas(fspace,in,deriv2);
        end
    end
    if size(G.P0,1)~=N || size(G.P0,2)~=N
        error('hjb:grid','Square collocation basis required');
    end
else
    G.P0=speye(N);
end
end
