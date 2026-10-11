function R=hjb_regimes(model)
if isfield(model,'nregimes') && ~isempty(model.nregimes)
    R=model.nregimes;
elseif isfield(model,'Q') && isnumeric(model.Q) && ~isempty(model.Q)
    if ndims(model.Q)==3, R=size(model.Q,2);
    else, R=size(model.Q,1); end
else
    R=1;
end
if ~isscalar(R) || R<1 || R~=floor(R)
    error('hjb:regimes','nregimes must be a positive integer');
end
end
