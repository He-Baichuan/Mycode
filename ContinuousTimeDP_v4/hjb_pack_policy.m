function x=hjb_pack_policy(p)
R=numel(p);
if R==1, x=p{1}; return; end
np=size(p{1},2);
if all(cellfun(@(z)size(z,2)==np,p))
    x=cat(3,p{:});
else
    x=p;
end
end
