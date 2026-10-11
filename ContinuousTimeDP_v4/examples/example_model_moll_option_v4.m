function model=example_model_moll_option_v4()
%EXAMPLE_MODEL_MOLL_OPTION_V4 Parameter match to Moll option_simple_LCP.m.
% Official script: rho=.05, x in [.1,1], flow sqrt(x), drift -.01,
% volatility .01*x, exercise reward 10, reflecting boundaries.
% Diffusion reflected FD stencil differs in boundary coefficient from
% Moll's one-sided corner convention; compare values after grid refinement.
model.nregimes=1; model.rho=.05;
model.boundary='reflecting';
model.stop=10;
model.candidates=@candidates;
end
function c=candidates(S,~,~,~)
N=size(S,1);
c=struct('x',zeros(N,1),'f',sqrt(S(:,1)),...
  'g',-.01*ones(N,1),'sigma',.01*S(:,1));
end
