function demo_skiba()
% Solve Moll/Kaplan Skiba butterfly-technology HJB using implicit FD.
[model,kmin,kmax,kstar]=example_model_skiba_moll();
F=fundefn('lin',230,kmin,kmax); k=funnode(F);
p=model.params{1}; V0=((k.^p.alpha).^(1-p.s)/(1-p.s))/p.rho;
opts=struct('maxit',350,'tol',1e-6,'res_tol',1e-5,...
    'Delta',1000,'showiters',true);
[V,c,res,info,A]=HJBsolve(model,F,k,V0,opts); %#ok<ASGLU>
figure;
subplot(2,1,1); plot(k,V,'LineWidth',1.5); hold on;
xline(kstar,'--'); title('Moll/Kaplan Skiba HJB: Value'); xlabel('Capital');
subplot(2,1,2); plot(k,c,'LineWidth',1.5); hold on;
xline(kstar,'--'); title('Optimal consumption'); xlabel('Capital');
fprintf('Skiba: converged %d residual %.3g\n',info.converged,info.residual_max);
end
