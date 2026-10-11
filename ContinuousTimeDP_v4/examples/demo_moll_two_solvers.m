function demo_moll_two_solvers()
% DEMO_MOLL_TWO_SOLVERS Two numerical methods on Moll baseline HJB.
% Call run_setup before this demo (or add root+examples directories).
close all;
model=example_model_moll();
I=200; amin=-0.02; amax=2; s=2; r=0.03; z=[0.1 0.2];
fFD=fundefn('lin',I,amin,amax); aFD=funnode(fFD);
V0FD=zeros(I,2);
for j=1:2
 V0FD(:,j)=((z(j)+r*aFD).^(1-s)/(1-s))/model.rho;
end
fdopts=struct('showiters',true,'maxit',180,'Delta',1000,...
 'tol',1e-6,'res_tol',1e-5);
[VFD,cFD,resFD,ifd,A]=HJBsolve(model,fFD,aFD,V0FD,fdopts);
fC=fundefn('spli',55,amin,amax);
aC=funnode(fC);
V0C=zeros(numel(aC),2);
for j=1:2
 V0C(:,j)=((z(j)+r*aC).^(1-s)/(1-s))/model.rho;
end
colopts=struct('showiters',true,'maxit',160,'tol',1e-6,...
  'res_tol',1e-5,'collocation_damping',0.7);
[cv,sC,VC,cC,resC,ic]=scsolve_jump(model,fC,aC,V0C,[],colopts);
figure;
subplot(1,2,1);
plot(aFD,VFD(:,1),aFD,VFD(:,2),'LineWidth',1.2);
xlabel('Assets'); ylabel('Value'); title('FD + Implicit');
subplot(1,2,2);
plot(sC,VC(:,1),sC,VC(:,2),'LineWidth',1.2);
xlabel('Assets'); ylabel('Value'); title('Collocation + PFI');
figure;
plot(aFD,cFD(:,1,1),aFD,cFD(:,1,2),... 
     sC,cC(:,1,1),'--',sC,cC(:,1,2),'--','LineWidth',1.1);
xlabel('Assets'); ylabel('Consumption');
legend('FD z_1','FD z_2','Colloc z_1','Colloc z_2','Location','best');
fprintf('FD: converged=%d residual=%.3g rowerr=%.3g\n',...
 ifd.converged,ifd.residual_max,max(abs(full(sum(A,2)))));
fprintf('Colloc: converged=%d residual=%.3g\n',ic.converged,ic.residual_max);
end
