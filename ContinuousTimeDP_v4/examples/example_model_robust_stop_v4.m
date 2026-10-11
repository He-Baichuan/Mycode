function model=example_model_robust_stop_v4()
%EXAMPLE_MODEL_ROBUST_STOP_V4 Two regimes, minimizer picks switching rate.
% Regime 1 stops at 15; regime 2 continues. V=(15,50/3).
model.nregimes=2; model.rho=.1;
model.boundary='state_constraint';
model.candidates=@candidates;
model.adversary=@adversary;
model.stop=[15;-inf];
end
function c=candidates(S,~,~,j)
N=size(S,1);
c=struct('x',zeros(N,1),'f',2*(j==2)*ones(N,1),...
 'g',zeros(N,1),'sigma',[]);
end
function nature=adversary(S,~,~,j,~)
N=size(S,1);
for t=1:2
    rate=.1*t;
    q=zeros(N,2);
    q(:,3-j)=rate;
    nature(t)=struct('q',q,'x',rate*ones(N,1)); %#ok<AGROW>
end
end
