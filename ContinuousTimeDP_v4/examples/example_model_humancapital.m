function [model,amin,amax,hmin,hmax]=example_model_humancapital()
%EXAMPLE_MODEL_HUMANCAPITAL Moll wealth/human capital model, HJB only.
% Source: https://benjaminmoll.com/codes/ human_capital.m
% Controls c (consumption), ell (education/learning time).
p=struct('gamma',2,'theta',.2,'r',.05,'rho',.05,'delta',.07,...
    'w',1,'alpha',.5);
ella=(p.alpha*p.delta)/(p.rho+p.delta);
hstar=(ella^p.alpha*p.theta/p.delta)^(1/(1-p.alpha));
amin=0; amax=5; hmin=.7*hstar; hmax=1.3*hstar;
p.lower=[amin hmin]; p.upper=[amax hmax];
model.nregimes=1; model.rho=p.rho;
model.params={p}; model.boundary={'state_constraint','state_constraint'};
model.candidates=@candidates;
end
function out=candidates(S,D,~,~,p)
a=S(:,1); h=S(:,2); N=numel(a);
% Optimize each upwind quadrant, including its zero-drift boundaries.
va={D.forward(:,1),D.backward(:,1)};
vh={D.forward(:,2),D.backward(:,2)};
ell0=((p.delta/p.theta).*h.^(1-p.alpha)).^(1/p.alpha);
out=repmat(struct('x',[],'f',[],'g',[],'sigma',[],'feasible',[]),1,4);
idx=0;
for ia=1:2
    for ih=1:2
        lo=zeros(N,1); hi=ones(N,1);
        if ih==1,lo=max(lo,ell0);else,hi=min(hi,ell0);end
        lowh=h<=p.lower(2)+1e-12; highh=h>=p.upper(2)-1e-12;
        lo(lowh)=max(lo(lowh),ell0(lowh));
        hi(highh)=min(hi(highh),ell0(highh));
        force=(ia==2 & a<=p.lower(1)+1e-12) | ...
            (ia==1 & a>=p.upper(1)-1e-12);
        objective=@(ell)hamiltonian(ell,a,h,va{ia},vh{ih},ia,force,p);
        [ell,H]=goldenx(objective,lo,hi);
        % Endpoint comparison also covers negative human-capital gradients.
        for endpoint={lo,hi}
            trial=endpoint{1}; score=objective(trial); take=score>H;
            ell(take)=trial(take); H(take)=score(take);
        end
        [~,c]=hamiltonian(ell,a,h,va{ia},vh{ih},ia,force,p);
        idx=idx+1; out(idx)=pack(c,ell,a,h,p);
    end
end
end
function [H,c]=hamiltonian(ell,a,h,va,vh,ia,force,p)
cash=p.r*a+p.w*h.*(1-ell);
c=max(va,1e-10).^(-1/p.gamma);
if ia==1,c=min(c,cash);else,c=max(c,cash);end
c(force)=cash(force); c=max(c,1e-12);
H=c.^(1-p.gamma)/(1-p.gamma)+va.*(cash-c)+...
    vh.*(p.theta*(ell.*h).^p.alpha-p.delta*h);
end
function q=pack(c,ell,a,h,p)
mu_a=p.r*a+p.w*h.*(1-ell)-c;
mu_h=p.theta*(ell.*h).^p.alpha-p.delta*h;
mu_a(abs(mu_a)<1e-12)=0; mu_h(abs(mu_h)<1e-12)=0;
valid=isfinite(c)&c>0&ell>=0&ell<=1;
q=struct('x',[c,ell],'f',c.^(1-p.gamma)/(1-p.gamma),...
    'g',[mu_a,mu_h],'sigma',[],'feasible',valid);
end
