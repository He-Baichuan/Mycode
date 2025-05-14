

var c k n y i g A r w Rk ;
varexo eA eg;

parameters alpha beta chi delta theta rho_A rho_G omega sigma_A sigma_G;

// 参数设定
alpha   = 1/3;
beta    = 0.99;
chi     = 1;
delta   = 0.025;
theta   = 4;
rho_A   = 0.95;
rho_G   = 0.97;
omega   = 0.2;
sigma_A = 1;
sigma_G = 1;
kN  = (alpha / (1/beta - (1 - delta)))^(1/(1 - alpha));
yN  = kN^alpha;
model;
    // 消费–资本 Euler 方程
    1/c = beta * (1/c(+1) * (Rk(+1) + 1 - delta));
    // 消费–资产 Euler 方程（可选，若有债券市场）
    1/c = beta * (1/c(+1) * (1 + r));

    // 边际产出价格
    w   = (1 - alpha) * A * k(-1)^alpha * n^(-alpha);
    Rk  = alpha * A * k(-1)^(alpha-1) * n^(1-alpha);

    // 劳动供给条件
    theta * n^chi = w / c;

    // 生产函数与资源约束
    y   = A * k(-1)^alpha * n^(1-alpha);
    y   = c + i + g;

    // 资本累积
    k = i + (1 - delta) * k(-1);

    // 政府支出与技术冲击过程
    log(g) = (1 - rho_G) * log(omega * y) + rho_G * log(g(-1)) + sigma_G * eg;
    log(A) = rho_A * log(A(-1))                     + sigma_A * eA;
end;

// 稳态初值计算
initval;
    A   = 1;
    // 资租率稳态
    Rk  = 1/beta - (1 - delta);
    
    // 工资率与产出
    w   = (1 - alpha) * kN^alpha;
    
    // 投资与消费
    i   = delta * kN;
    c   = (1 - omega) * yN + delta * kN;
    // 劳动供给
    n   = ((1/theta) * w / c)^(1/(1 + chi));
    // 绝对水平
    k   = kN * n;
    y   = yN * n;
    g   = omega * y;
end;

shocks;
    var eg; stderr sigma_G^2;
    var eA; stderr sigma_A^2;
end;

steady;
stoch_simul(order=1,irf = 180);
