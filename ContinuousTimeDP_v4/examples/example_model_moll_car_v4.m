function model=example_model_moll_car_v4()
%EXAMPLE_MODEL_MOLL_CAR_V4 HJB/QVI adapter for Moll car.m (not GE).
% The wealth grid needs amin=-.02, amax=3. Regime 1: no car,
% regime 2: car. Buying pays .2 assets, selling receives .1 assets,
% flow benefit car=+.25. Income y=.1, r=.045, rho=.05, CRRA s=2.
% This is a model adapter, not yet MATLAB-validated against Moll output.
model.nregimes=2; model.rho=.05;
model.boundary='state_constraint';
model.params={2,.045,.1,.25,.2,.1,-.02,3};
model.candidates=@candidates;
model.intervention=@intervention;
end
function cs=candidates(S,D,~,j,s,r,y,kappa,~,~,amin,amax)
N=size(S,1); cash=y+r*S(:,1);
C=[D.forward(:,1),D.backward(:,1),cash.^(-s)];
cs=repmat(struct('x',[],'f',[],'g',[],'sigma',[]),1,3);
for t=1:3
    cp=max(C(:,t),1e-11).^(-1/s);
    % No outward savings at asset state constraints.
    atLower=S(:,1)<=amin+1e-10;
    atUpper=S(:,1)>=amax-1e-10;
    cp(atLower)=min(cp(atLower),cash(atLower));
    cp(atUpper)=max(cp(atUpper),cash(atUpper));
    cp=max(cp,1e-10);
    cs(t).x=cp;
    cs(t).f=cp.^(1-s)/(1-s) + kappa*(j==2);
    cs(t).g=cash-cp;
    cs(t).sigma=[];
end
end
function op=intervention(S,j,~,~,s,r,y,kappa,p0,p1,amin,amax) %#ok<INUSD>
N=size(S,1); a=S(:,1);
if j==1
    % Infeasible buyers get a benign in-domain target but are rejected.
    target=max(amin,a-p0);
    feasible=(a-p0)>=amin;
    op=struct('target',target,'to',2,'reward',0,...
              'feasible',feasible,'action',ones(N,1));
else
    % Moll caps the sale receipt if it would exceed amax.
    target=min(amax,a+p1);
    op=struct('target',target,'to',1,'reward',0,...
              'feasible',true(N,1),'action',2*ones(N,1));
end
end
