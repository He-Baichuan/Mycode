function res = rbc_resid(m, X, E)
     % 从 m.par 里取出所有参数
  alpha = m.par.alpha;
  beta  = m.par.beta;
  gamma = m.par.gamma; 
  delta = m.par.delta;
  rho   = m.par.rho;

 T = m.T;
 p = m.nx;
  res = zeros(T,p);

    Xmat = reshape(X, m.nx, m.T).';
    Z = Xmat(:, m.idx_now('Z'));
    R = Xmat(:, m.idx_now('R'));
    K = Xmat(:, m.idx_now('K'));
    Y = Xmat(:, m.idx_now('Y'));
    C = Xmat(:, m.idx_now('C'));
    Z_l = ModelUtils.lag(Z);
    K_l = ModelUtils.lag(K);
    R_p= ModelUtils.lead(R);
    C_p = ModelUtils.lead(C);



    eps = E;
    
        res(:,1) = -C.^(-gamma) + beta.*R_p.*C_p.^(-gamma);
        res(:,2) = -R + alpha*Z.*K_l.^(alpha-1) + 1 - delta;
        res(:,3) = -K+ (1-delta)*K_l + Y - C;
        res(:,4) = -Y + Z.*K_l.^alpha;
        res(:,5) = -log(Z) + rho*log(Z_l) + eps;
        res = res';
        res = res(:);


end