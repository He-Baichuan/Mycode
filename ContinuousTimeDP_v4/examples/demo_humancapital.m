function demo_humancapital()
% Two-dimensional controlled HJB with constrained upwind controls.
[model,amin,amax,hmin,hmax]=example_model_humancapital();
F=fundefn('lin',[36 23],[amin hmin],[amax hmax]);
nodes=funnode(F); G=hjb_grid(F,nodes,'fd');
p=model.params{1}; a=G.S(:,1); h=G.S(:,2);
ell0=min(1,max(0,((p.delta/p.theta).*h.^(1-p.alpha)).^(1/p.alpha)));
c0=p.r*a+p.w*h.*(1-ell0);
V0=c0.^(1-p.gamma)/(1-p.gamma)/p.rho;
opts=struct('maxit',300,'Delta',100,'showiters',true,...
  'tol',1e-6,'res_tol',2e-5);
[V,pol,res,info]=HJBsolve(model,F,nodes,V0,opts); %#ok<ASGLU>
figure;
subplot(1,2,1); imagesc(nodes{1},nodes{2},reshape(V,G.n)'); axis xy;
xlabel('Assets');ylabel('Human capital');title('Value');colorbar;
subplot(1,2,2); imagesc(nodes{1},nodes{2},reshape(pol(:,2),G.n)'); axis xy;
xlabel('Assets');ylabel('Human capital');title('Education time');colorbar;
fprintf('Human capital: converged=%d residual=%.3g\n',...
    info.converged,info.residual_max);
end
