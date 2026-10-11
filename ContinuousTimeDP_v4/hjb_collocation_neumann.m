function [M,b]=hjb_collocation_neumann(M,b,G,R,model)
% Apply the same reflecting rows in policy evaluation and residual checks.
if ~isfield(model,'boundary') || isempty(model.boundary), return; end
bc=model.boundary;
if ischar(bc), bc=repmat({bc},1,G.d); end
replaced=false(size(b,1),1);
for k=1:G.d
    if ~strcmpi(bc{k},'reflecting'), continue; end
    nodes=find(G.subs(:,k)==1|G.subs(:,k)==G.n(k));
    for reg=1:R
        rr=(reg-1)*G.N+nodes;
        valid=~replaced(rr); rr=rr(valid); ij=nodes(valid);
        if isempty(rr), continue; end
        M(rr,:)=0;
        cols=(reg-1)*G.N+(1:G.N);
        M(rr,cols)=G.P1{k}(ij,:);
        b(rr)=0;
        replaced(rr)=true;
    end
end
end
