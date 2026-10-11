function [model,kmin,kmax,kstar]=example_model_skiba_moll()
%EXAMPLE_MODEL_SKIBA_MOLL Production specification from Moll/Kaplan Skiba demo.
% SOURCE: https://benjaminmoll.com/codes/ HJB_NGM_skiba.m
p=struct('s',2,'alpha',.3,'rho',.05,'delta',.05,...
    'AH',.6,'AL',.4,'kappa',2);
kssH=(p.alpha*p.AH/(p.rho+p.delta))^(1/(1-p.alpha))+p.kappa;
kstar=p.kappa/(1-(p.AL/p.AH)^(1/p.alpha));
kmin=.001*kssH; kmax=1.3*kssH;
model.nregimes=1; model.rho=p.rho;
model.boundary={'state_constraint'};
model.params={p}; model.candidates=@skiba_candidates;
end
function cand=skiba_candidates(S,D,~,~,p)
k=S(:,1); N=numel(k);
y=max(p.AH*max(k-p.kappa,0).^p.alpha,p.AL*k.^p.alpha);
cash=y-p.delta*k;
if any(cash<=0), error('Skiba grid must have positive net output'); end
gradients=[D.forward(:,1),D.backward(:,1),max(cash,1e-8).^(-p.s)];
cand=repmat(struct('x',[],'f',[],'g',[],'sigma',[]),1,3);
for idx=1:3
    c=max(gradients(:,idx),1e-10).^(-1/p.s);
    if idx==3, c=cash; end
    c(1)=min(c(1),cash(1)); c(end)=max(c(end),cash(end));
    cand(idx).x=c;
    cand(idx).f=c.^(1-p.s)/(1-p.s);
    cand(idx).g=cash-c;
    cand(idx).sigma=[];
end
end
