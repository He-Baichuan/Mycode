function model=example_model_optstop_v4()
%EXAMPLE_MODEL_OPTSTOP_V4 Deterministic obstacle with exact kink at a=0.5.
% rho=0.1, f(a)=2*a, zero drift: V=max(10,20*a).
model.nregimes=1; model.rho=.1;
model.boundary='state_constraint';
model.candidates=@candidates;
model.stop=10;
end
function c=candidates(S,~,~,~)
N=size(S,1);
c=struct('x',zeros(N,1),'f',2*S(:,1),'g',zeros(N,1),'sigma',[]);
end
