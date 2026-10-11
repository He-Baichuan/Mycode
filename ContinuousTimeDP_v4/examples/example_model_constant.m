function model=example_model_constant(varargin)
%EXAMPLE_MODEL_CONSTANT Analytic benchmark for Markov switching and drift.
% Each regime receives spatially constant flow reward f_j; solution
% V_j = [(rho*I-Q)\f]_j is constant in continuous states.
R=2;
if nargin>=1, R=varargin{1}; end
model.nregimes=R; model.rho=0.07;
if R==2, model.Q=[-0.1 0.1;0.04 -0.04];
else, model.Q=zeros(R,R); end
model.boundary={'state_constraint'};
model.candidates=@constant_candidates;
end
function c=constant_candidates(S,~,~,j)
N=size(S,1); d=size(S,2);
g=zeros(N,d);
for k=1:d, g(:,k)=0.2*(0.5-S(:,k)); end
c=struct('x',zeros(N,1),'f',ones(N,1)*j,'g',g,'sigma',[]);
end
