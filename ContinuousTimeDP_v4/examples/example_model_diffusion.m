function model=example_model_diffusion()
%EXAMPLE_MODEL_DIFFUSION Reflecting diffusions; exact V=flow/rho.
model.nregimes=1; model.rho=0.1;
model.boundary={'reflecting'};
model.candidates=@candidate;
end
function c=candidate(S,~,~,~)
N=size(S,1);
c=struct('x',zeros(N,1),'f',2*ones(N,1),...
    'g',0.1*(0.5-S),'sigma',sqrt(0.04)*ones(N,1));
end
