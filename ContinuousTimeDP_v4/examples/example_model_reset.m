function model=example_model_reset()
%EXAMPLE_MODEL_RESET Poisson reset to a=.5 at rate .4; exact affine V.
model.nregimes=1; model.rho=0.1;
model.boundary={'state_constraint'};
model.candidates=@reset_candidates;
model.jump=@reset_jump;
end
function c=reset_candidates(S,~,~,~)
N=size(S,1);
c=struct('x',zeros(N,1),'f',S(:,1),'g',zeros(N,1),'sigma',[]);
end
function channels=reset_jump(S,~,~)
channels=struct('rate',0.4,'target',0.5*ones(size(S)),...
    'to',1,'weight',1,'reward',0);
end
