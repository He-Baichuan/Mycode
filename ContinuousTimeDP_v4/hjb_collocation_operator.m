function L=hjb_collocation_operator(G,mu,sigma)
%HJB_COLLOCATION_OPERATOR Map basis coefficients to generator on nodes.
N=G.N; d=G.d; L=sparse(N,N);
for k=1:d
    L=L+spdiags(mu(:,k),0,N,N)*G.P1{k};
end
if isempty(sigma), return; end
if d==1 && isvector(sigma), sigma=reshape(sigma,N,1,1); end
if ismatrix(sigma) && d>1, sigma=reshape(sigma,N,d,1); end
q=size(sigma,3);
for k=1:d
    for l=k:d
        aij=sum(reshape(sigma(:,k,:),N,q).*reshape(sigma(:,l,:),N,q),2);
        if k==l, aij=aij/2; end
        if any(aij~=0)
            L=L+spdiags(aij,0,N,N)*G.P2{k,l};
        end
    end
end
end
