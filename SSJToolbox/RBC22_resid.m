function res = RBC22_resid(m, X, E)
% 从 m.par 里取出所有参数
  alpha = m.par.alpha;
  beta  = m.par.beta;
  theta = m.par.theta;
  delta = m.par.delta;
  rho   = m.par.rho;

 T = m.T;
 p = m.nx;
  res = zeros(p*T,1);
  for t=1:T
    idx0 = p*(t-1);
    % 当前、滞后、前导期逻辑同前……
    a = X(idx0+1);  
    r = X(idx0+2);
    k = X(idx0+3);  
    y = X(idx0+4);
    c = X(idx0+5);
    Rk = X(idx0+6);
    w = X(idx0+7);  
    i = X(idx0+8);
    n = X(idx0+9);
    if t==1
      a_l = a;  
      k_l = k;
    else
      a_l = X(idx0-p+1);
      k_l = X(idx0-p+3);
    end
    if t==T
      Rk_p = Rk;  
      c_p = c;
      % a_p = a;
      % k_p = k;
    else
      Rk_p = X(idx0+p+6);
      c_p = X(idx0+p+5);
      % a_p = X(idx0+p+1);
      % k_p = X(idx0+p+3);
    end
    eps = E(t);
    % 写五条方程
    res(idx0+1) = -(1/c-beta*((1/c_p)*(Rk_p+1-delta)));
    res(idx0+2) = -(theta/(1-n)-w/c);
    res(idx0+3) = -(k-a*k_l^(alpha)*n^(1-alpha)+c-(1-delta)*k_l);
    res(idx0+4) = -(log(a)-rho*log(a_l)-eps);
    res(idx0+5) = -(y-c-i);
    res(idx0+6) = -(y-a*k_l^(alpha)*n^(1-alpha));
    res(idx0+7) = -(1/c-beta*((1/c_p)*(1+r)));
    res(idx0+8) = -(w-(1-alpha)*a*k_l^(alpha)*n^(-alpha));
    res(idx0+9) = -(Rk-alpha*a*k_l^(alpha-1)*n^(1-alpha));
  end

end