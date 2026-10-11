function demo_moll_D_v4()
%DEMO_MOLL_D_V4 Experimental (not verified) Moll D1/D3 HJB adapters.
% D1: option_simple_LCP.m — FD reflecting obstacle.
% D3: car.m — FD 2-regime QVI with wealth reset.
run_setup();
opts=struct('maxit',200,'tol',1e-7,'res_tol',1e-5,...
    'showiters',true,'Delta',1000,'lcp_tol',1e-6);
% Option model
m=example_model_moll_option_v4();
fspace=fundefn('lin',151,.1,1); x=funnode(fspace);
[V,~,R,info,A]=HJBsolve(m,fspace,x,10*ones(151,1),opts); %#ok<ASGLU>
fprintf('Moll option: converged=%d, residual=%.3e, stop nodes=%d\n',...
    info.converged,info.residual_max,nnz(info.stop_mask));
figure(1); plot(x,V,x,10*ones(size(x)),'--'); grid on;
title('Moll option HJB/VI adapter'); legend('value','exercise reward');
% Car model
m=example_model_moll_car_v4();
fspace=fundefn('lin',201,-.02,3); a=funnode(fspace);
V0=repmat((-1./(.1+.045*a))/.05,1,2);
[V,~,R,info,A]=HJBsolve(m,fspace,a,V0,opts); %#ok<ASGLU>
fprintf('Moll car: converged=%d, residual=%.3e, buy nodes=%d, sell nodes=%d\n',...
    info.converged,info.residual_max,...
    nnz(info.intervention_mask(:,1)),nnz(info.intervention_mask(:,2)));
figure(2); plot(a,V); grid on;
title('Moll car HJB/QVI adapter'); legend('without car','with car');
end
