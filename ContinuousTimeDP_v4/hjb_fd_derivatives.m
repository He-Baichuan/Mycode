function D=hjb_fd_derivatives(G,W)
%HJB_FD_DERIVATIVES One-sided derivatives for possibly NONUNIFORM grids.
N=G.N; d=G.d;
if numel(W)~=N, error('hjb:derivatives','Value dimension mismatch'); end
D.forward=zeros(N,d); D.backward=zeros(N,d);
for k=1:d
    idx=G.subs(:,k); validF=idx<G.n(k); validB=idx>1;
    nxt=find(validF); prv=find(validB);
    hf=diff(G.nodes{k});
    D.forward(nxt,k)=(W(nxt+G.stride(k))-W(nxt))./hf(idx(nxt));
    D.backward(prv,k)=(W(prv)-W(prv-G.stride(k)))./hf(idx(prv)-1);
    upper=find(~validF); lower=find(~validB);
    D.forward(upper,k)=D.backward(upper,k);
    D.backward(lower,k)=D.forward(lower,k);
end
D.central=0.5*(D.forward+D.backward);
D.grad=D.central;
end
