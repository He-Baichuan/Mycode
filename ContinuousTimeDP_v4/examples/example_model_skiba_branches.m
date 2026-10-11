function model=example_model_skiba_branches()
%EXAMPLE_MODEL_SKIBA_BRANCHES Stylized nonconvex technology switching.
% A benchmark for GLOBAL branch choice; not Moll's exact NGM Skiba model.
% f1(a)=1+.8*a ; f2(a)=1.8-.8*a, g=0; kink in V at a=0.5.
model.nregimes=1; model.rho=0.1;
model.boundary={'state_constraint'};
model.candidates=@candidates;
end
function cs=candidates(S,~,~,~)
N=size(S,1);
cs(1)=struct('x',ones(N,1),'f',1+0.8*S(:,1),...
  'g',zeros(N,1),'sigma',[]);
cs(2)=struct('x',2*ones(N,1),'f',1.8-0.8*S(:,1),...
  'g',zeros(N,1),'sigma',[]);
end
