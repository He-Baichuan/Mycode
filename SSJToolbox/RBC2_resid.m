function res = RBC2_resid(m, X, E)
% 从 m.par 里取出所有参数
 %  alpha = m.par.alpha;
 %  beta  = m.par.beta;
 %  theta = m.par.theta;
 %  delta = m.par.delta;
 %  rho   = m.par.rho;
 % 
 T = m.T;
 p = m.nx;
  res = zeros(T,p);
 %  % for t=1:T
 %  %   idx0 = p*(t-1);
 %  %   % 当前、滞后、前导期逻辑同前……

 %  a = ModelUtils.contemp(X,'a',m);
 %  r = ModelUtils.contemp(X,'r',m);
 %  k = ModelUtils.contemp(X,'k',m);
 %  y = ModelUtils.contemp(X,'y',m);
 %  c = ModelUtils.contemp(X,'c',m);
 %  Rk = ModelUtils.contemp(X,'Rk',m);
 %  w = ModelUtils.contemp(X,'w',m);
 %  i = ModelUtils.contemp(X,'i',m);
 %  n = ModelUtils.contemp(X,'n',m);
 %  %   if t==1
 %  %     a_l = a;  
 %  %     k_l = k;
 %  %   else
 %  %     a_l = X(idx0-p+1);
 %  %     k_l = X(idx0-p+3);
 %  %   end
 %  a_l = ModelUtils.lag(a);
 %  k_l = ModelUtils.lag(k);
 %  %   if t==T
 %  %     Rk_p = Rk;  
 %  %     c_p = c;
 %  %     % a_p = a;
 %  %     % k_p = k;
 %  %   else
 %  %     Rk_p = X(idx0+p+6);
 %  %     c_p = X(idx0+p+5);
 %  %     % a_p = X(idx0+p+1);
 %  %     % k_p = X(idx0+p+3);
 %  %   end
 %  Rk_p = ModelUtils.lead(Rk,T);
 %  c_p = ModelUtils.lead(c,T);
 % 
 % 
    eps = E;
 %    % 写五条方程
    
 %  % end
 %    res = res';
 %    res = res(:);
% 矢量化残差：一次性返回 (T·n_eq)×1 列向量
% -------------------------------------------------------------
    % 快速拆分 X：reshape 后列即时间
    Xmat = reshape(X, m.nx, m.T).';         % T×nx
    a = Xmat(:, m.idx_now('a'));
    r = Xmat(:, m.idx_now('r'));
    k = Xmat(:, m.idx_now('k'));
    y = Xmat(:, m.idx_now('y'));
    c = Xmat(:, m.idx_now('c'));
    Rk= Xmat(:, m.idx_now('Rk'));
    w = Xmat(:, m.idx_now('w'));
    i = Xmat(:, m.idx_now('i'));
    n = Xmat(:, m.idx_now('n'));

    % lag / lead（用稳态补端点，见前文）
    a_l = ModelUtils.lag(a);
    k_l = ModelUtils.lag(k);
    Rk_p= ModelUtils.lead(Rk);
    c_p = ModelUtils.lead(c);

    % -------------- 计算九条方程 ----------------
    alpha = m.par.alpha;  beta = m.par.beta; theta = m.par.theta;
    delta = m.par.delta;  rho  = m.par.rho;

    
    res(:,1) = -(1./c-beta*((1./c_p).*(Rk_p+1-delta)));
    res(:,2) = -(theta./(1-n)-w./c);
    res(:,3) = -(k-a.*k_l.^(alpha).*n.^(1-alpha)+c-(1-delta).*k_l);
    res(:,4) = -(log(a)-rho*log(a_l)-eps);
    res(:,5) = -(y-c-i);
    res(:,6) = -(y-a.*k_l.^(alpha).*n.^(1-alpha));
    res(:,7) = -(1./c-beta*((1./c_p).*(1+r)));
    res(:,8) = -(w-(1-alpha)*a.*k_l.^(alpha).*n.^(-alpha));
    res(:,9) = -(Rk-alpha*a.*k_l.^(alpha-1).*n.^(1-alpha));
    res = res';                         % → 9×T
    res = res(:);                       % → (9T)×1
end


