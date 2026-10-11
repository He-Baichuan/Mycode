function model=example_model_qvi_v4()
%EXAMPLE_MODEL_QVI_V4 Nonlocal impulse regime switch (stylized car purchase).
% Regime 1 no flow, regime 2 flow 2, rho=.1, no drift.
% Regime 1 can pay 3 to move to regime 2, reset a to 1-a.
% Solution V_1(a)=17, V_2(a)=20 everywhere.
model.nregimes=2; model.rho=.1;
model.boundary='state_constraint';
model.candidates=@candidates;
model.intervention=@intervene;
end
function c=candidates(S,~,~,j)
N=size(S,1);
c=struct('x',zeros(N,1),'f',2*(j==2)*ones(N,1),...
    'g',zeros(N,size(S,2)),'sigma',[]);
end
function op=intervene(S,j,~,evalfun)
% evalfun allows callback authors to optimize target continuously if desired.
if j==2, op=[]; return; end
% demonstrate value interpolation to a nonlocal reset location
T=1-S;
if ~isequal(size(evalfun(T,2)),[size(S,1),1])
    error('hjb:test','Bad intervention evaluator');
end
op=struct('target',T,'to',2,'reward',-3,'action',ones(size(S,1),1));
end
