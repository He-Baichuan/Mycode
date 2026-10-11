function model=example_model_moll()
%EXAMPLE_MODEL_MOLL CRRA agent, asset state and two Poisson income regimes.
% Moll's s=2,r=.03,rho=.05,z=(.1,.2),la=(.02,.03).
model.nregimes=2; model.rho=0.05;
model.Q=[-0.02,0.02;0.03,-0.03];
model.boundary={'state_constraint'};
model.params={2,0.03,[.1 .2]};
model.candidates=@moll_candidates;
end

function cs=moll_candidates(S,D,~,reg,s,r,z)
N=size(S,1); cash=z(reg)+r*S(:,1);
grads=[D.forward(:,1),D.backward(:,1),D.central(:,1)];
cs=repmat(struct('x',[],'f',[],'g',[],'sigma',[]),1,3);
for k=1:3
    % CRRA FOC: c = (V_a)^(-1/s). Clamp only to prevent numerical
    % evaluation failure; boundary feasibility is enforced by solver.
    dv=max(grads(:,k),1e-12);
    c=min(max(dv.^(-1/s),1e-10),1e5);
    % At exact state constraints: prohibit asset drift out of support.
    atLower=S(:,1)<=min(S(:,1))+1e-12;
    atUpper=S(:,1)>=max(S(:,1))-1e-12;
    c(atLower)=min(c(atLower),cash(atLower));
    c(atUpper)=max(c(atUpper),cash(atUpper));
    c=max(c,1e-10);
    cs(k).x=c;
    cs(k).f=(c.^(1-s))/(1-s);
    cs(k).g=cash-c;
    cs(k).sigma=[];
end
end
