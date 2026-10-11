function demo_v4()
%DEMO_V4 Run v4 steady HJB-VI/QVI on finite difference and collocation.
run_setup();
mods={@example_model_optstop_v4,@example_model_qvi_v4,@example_model_robust_stop_v4};
labels={'Optimal stopping','Nonlocal impulse QVI','Robust optimal stopping'};
fd=struct('maxit',100,'tol',1e-8,'res_tol',1e-7,'showiters',false,...
    'Delta',1000,'lcp_backend','compecon','lcp_tol',1e-7);
co=struct('maxit',50,'tol',1e-8,'res_tol',1e-7,'showiters',false,...
    'lcp_backend','compecon','lcp_tol',1e-7);
for t=1:numel(mods)
    model=mods{t}();
    F=fundefn('lin',81,0,1); a=funnode(F);
    [V,~,~,id]=HJBsolve(model,F,a,[],fd);
    FC=fundefn('cheb',31,0,1); ac=funnode(FC);
    [~,~,W,~,~,ic]=scsolve_jump(model,FC,ac,[],[],co);
    fprintf('%s: FD residual %.2e, PFI residual %.2e, FD converged %d, Coll converged %d\n',...
        labels{t},id.residual_max,ic.residual_max,id.converged,ic.converged);
    figure(t); plot(a,V,'LineWidth',1.8); hold on;
    plot(ac,W,'--','LineWidth',1.5); hold off; grid on;
    title(labels{t}); xlabel('continuous state'); ylabel('value');
end
end
