function model=example_model_branches()
%EXAMPLE_MODEL_BRANCHES Two discrete choices; known optimal reward max.
model.nregimes=1; model.rho=0.1;
model.boundary={'state_constraint'};
model.candidates=@branches;
end
function out=branches(S,~,~,~)
N=size(S,1); d=size(S,2);
out(1)=struct('x',ones(N,1),'f',ones(N,1),...
    'g',zeros(N,d),'sigma',[]);
out(2)=struct('x',2*ones(N,1),'f',2*ones(N,1),...
    'g',zeros(N,d),'sigma',[]);
end
