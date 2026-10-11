# ContinuousTimeDP v4 `HJBsolve` 使用文档

本文说明 `D:/ContinuousTimeDP_v4/ContinuousTimeDP_v4/HJBsolve.m` 的使用方法。
`HJBsolve` 是一个稳态、无限期连续时间 HJB 有限差分求解器，支持：

- 普通 HJB：`max_u { f + A^u V }`；
- 稳健 HJB：`max_u min_h { f + A^(u,h) V }`；
- 最优停止变分不等式（VI）；
- 带内生干预的脉冲控制拟变分不等式（QVI）。

它只求解个体的 HJB/VI/QVI，不求解 KFE、稳态分布、一般均衡或转移动态。

## 1. 安装和路径

需要 MATLAB R2023b 或兼容版本，以及 Miranda--Fackler CompEcon 的 `CEtools`。
至少需要 `fundefn`、`funnode`、`funbas`、`gridmake`、`lcpsolve`；人力资本示例还使用 `goldenx`。

```matlab
addpath(genpath('D:/code/Collocation/CEtools'));
addpath('D:/ContinuousTimeDP_v4/ContinuousTimeDP_v4');
run_setup;

which HJBsolve -all
which hjb_grid -all
which lcpsolve -all
```

不要把 v3 和 v4 同时放在 MATLAB 路径中而不检查优先级。`which HJBsolve -all` 应该把 v4 目录列为第一项。

## 2. 调用格式

```matlab
[V, x, resid, info, A] = HJBsolve(model, fspace, snodes, V0, options);
```

参数含义：

| 参数 | 含义 |
|---|---|
| `model` | 模型结构体，至少包含 `rho`、`boundary`、`candidates` |
| `fspace` | CompEcon 的状态空间；FD 必须使用 `lin` 基函数 |
| `snodes` | 状态节点。传 `[]` 时自动使用 `funnode(fspace)`；多维时可传 cell 数组 |
| `V0` | 初始价值函数，传 `[]` 使用零初值；形状为 `N`×`R` 或等长向量 |
| `options` | 数值选项结构体，传 `[]` 使用默认值 |

其中 `N` 是连续状态节点数，`R` 是离散 regime 数。多维状态的节点顺序遵循 `gridmake`：第一维变化最快。
每个维度至少需要 3 个严格递增节点。

## 3. 最小可运行例子

v4 自带一个有解析解的最优停止例子：

```matlab
run_setup;

F = fundefn('lin', 61, 0, 1);
a = funnode(F);
model = example_model_optstop_v4();

opts = struct('maxit',120, 'Delta',1000, ...
    'tol',1e-9, 'res_tol',1e-7, 'showiters',true);

[V, x, resid, info, A] = HJBsolve(model, F, a, [], opts);

fprintf('converged=%d, iterations=%d, residual=%.3e\n', ...
    info.converged, info.iterations, info.residual_max);
plot(a,V,'LineWidth',1.4); hold on;
plot(a,10*ones(size(a)),'--'); grid on;
legend('V','stopping payoff','Location','best');
```

该例子的模型是 `rho=.1`、流收益 `f(a)=2a`、停止收益 `G=10`，解析解为
`V(a)=max(10,20a)`，停止边界为 `a=.5`。

## 4. `model` 结构体

### 4.1 必需字段

```matlab
model.nregimes = 1;                 % 可省略，默认 1
model.rho = .05;                    % 标量，或由 hjb_eval_rho 扩展的模型设置
model.boundary = 'state_constraint';
model.candidates = @candidate_callback;
```

`boundary` 可以是一个字符串，也可以按连续状态维度给出 cell：

```matlab
model.boundary = {'state_constraint','reflecting'};
```

`state_constraint` 表示状态约束边界：控制不能产生向计算域外的漂移，法向扩散必须为零。
`reflecting` 表示反射边界，FD 算子使用镜像 ghost-node 的 Neumann 处理。

### 4.2 控制候选回调

求解器在每个节点、每个 regime 调用：

```matlab
candidates = model.candidates(S, D, W, j, params{:});
```

其中：

- `S`：`N`×`d` 的状态节点；
- `D.forward`、`D.backward`：各维的一阶前向、后向差分；
- `D.central`、`D.grad`：中心差分；
- `D.value`：当前 regime 的节点价值；
- `D.all_values`：所有 regime 的节点价值；
- `W`：`N`×`R` 的当前节点价值；
- `j`：当前 regime 编号；
- `params`：`model.params` 中的参数 cell 内容。

回调必须返回一个或多个结构体候选。每个候选至少包含：

```matlab
candidate.x       % N-by-p，控制变量
candidate.f       % N-by-1，流收益
candidate.g       % N-by-d，连续状态漂移
```

可选字段：

```matlab
candidate.sigma   % N-by-d-by-q，Brownian loading
candidate.feasible% N-by-1 logical，可行性掩码
candidate.q       % regime 转移率
candidate.jump    % 有限强度跳跃通道
```

返回多个候选时，求解器在每个节点比较完整的离散 Hamiltonian，选择最大者；因此非凸控制或多个局部分支应全部返回，不能只返回一个局部一阶条件解。

一个最简单的单 regime 模型：

```matlab
model.nregimes = 1;
model.rho = .1;
model.boundary = 'state_constraint';
model.params = {};
model.candidates = @constant_candidate;

function c = constant_candidate(S,~,~,~)
N = size(S,1);
c = struct('x',zeros(N,1), ...
    'f',ones(N,1), ...
    'g',zeros(N,1), ...
    'sigma',[]);
end
```

### 4.3 Markov regime switching

固定转移矩阵可以直接写入 `model.Q`：

```matlab
model.nregimes = 2;
model.Q = [-.02 .02; .03 -.03];
```

对角线会被求解器置为当前 regime 的离开率约定；非对角元素必须非负。
也可以在候选中提供 `candidate.q`，其形状可以是 `R`×`R`、`N`×`R`，或 `N`×`R`×`R` 的逐节点形式。

### 4.4 稳健控制

稳健模型增加：

```matlab
model.adversary = @nature_callback;
```

回调格式：

```matlab
variants = model.adversary(S,D,W,j,candidate,params{:});
```

它返回多个 Nature 候选。每个候选可以覆盖控制候选的 `f`、`g`、`sigma`、`q`、`jump` 或 `feasible`，并可用 `.x` 保存 Nature 的扭曲变量。求解器计算的是

```text
max over controller candidates, then min over supplied Nature candidates
```

这个最小化只对回调提供的候选成立。连续的 Nature 优化必须在回调内部解析或数值求解。

### 4.5 最优停止

设置停止收益：

```matlab
model.stop = 10;                 % 所有节点、所有 regime
```

也可以使用 regime 向量、`N`×`R` 数组或函数：

```matlab
model.stop = [15; -inf];
% 或
model.stop = @(S,j,p) 15*(j==1)*ones(size(S,1),1);
```

`-Inf` 表示该 regime 或节点禁止停止。求解器使用 LCP 处理

```text
V >= G,
F = (rho + 1/Delta) V - A V - f - V_old/Delta,
V = G => F >= 0,
V > G => F = 0.
```

### 4.6 脉冲控制 QVI

设置：

```matlab
model.intervention = @intervention_callback;
```

回调格式：

```matlab
options = model.intervention(S,j,W,evalfun,params{:});
```

每个干预候选可包含：

```matlab
option.target    % N-by-d，干预后的连续状态
option.to        % 目标 regime，标量或 N-vector
option.reward    % 即时净收益；成本应写成负数
option.feasible  % N-vector logical
option.action    % N-by-p，保存干预控制
```

`evalfun(T,to)` 会在当前价值函数上评价任意目标状态。求解器形成

```text
M V(s,j) = max_xi { reward(s,j,xi) + V(target(s,j,xi),to) }
```

并求解 `max(MV-V, HJB)=0`。通用求解器不会自动在无限连续干预集合中寻找全局最优目标；模型回调必须提供足够的候选，或在回调内部完成目标优化。

## 5. 有限差分网格和扩散

一维：

```matlab
F = fundefn('lin', 501, 0, 10);
nodes = funnode(F);
```

二维：

```matlab
F = fundefn('lin',[101 61],[0 .5],[5 2]);
nodes = funnode(F);       % cell 数组
G = hjb_grid(F,nodes,'fd');
```

非均匀网格也可以使用。若有非零交叉扩散，当前正型差分模板要求各维使用均匀网格，并且协方差满足正型条件；否则会报 `hjb:cross` 或 `hjb:covariance`。

扩散 loading 的格式是 `N`×`d`×`q`。例如单一收入 Brownian：

```matlab
sigma = zeros(N,2,1);
sigma(:,2,1) = .2*S(:,2);
```

法向扩散存在时，不应把对应边界设置成 `state_constraint`；应使用 `reflecting` 或重新定义模型边界。

## 6. 数值选项

常用选项及默认值：

| 选项 | 默认值 | 作用 |
|---|---:|---|
| `maxit` | 200 | 最大政策迭代次数 |
| `Delta` | 1000 | 隐式伪时间步长；越大越接近直接稳态更新，过大可能影响数值稳定性 |
| `tol` | 1e-7 | 相邻价值函数的最大变化容差 |
| `res_tol` | 1e-6 | 更新后 HJB/VI/QVI 残差容差 |
| `stop_on_residual` | true | 是否要求同时满足 `tol` 和 `res_tol` |
| `showiters` | true | 是否打印每次迭代诊断 |
| `check_generator` | true | 是否检查生成矩阵行和与非对角元正性 |
| `generator_tol` | 1e-9 | 生成矩阵检查容差 |
| `infeasible_tol` | 1e-11 | 边界漂移和扩散的可行性容差 |
| `qvi_damping` | 1 | QVI 更新阻尼，取值 `(0,1]` |
| `lcp_backend` | `'compecon'` | LCP 后端，可选 `'compecon'` 或 `'moll'` |
| `lcp_tol` | 1e-7 | LCP 互补条件容差 |
| `obstacle_tol` | 1e-6 | 判断活动障碍的容差 |
| `jump_outside` | `'error'` | 跳跃目标出界时的处理，可选 `'clip'` |

典型配置：

```matlab
opts = struct('maxit',300, 'Delta',100, ...
    'tol',1e-7, 'res_tol',1e-6, ...
    'showiters',true, 'qvi_damping',.7);
```

`info.converged=true` 的含义是：最后一次迭代同时满足价值函数变化和残差要求；它不等于误差有理论上界，也不替代网格加密检查。

## 7. 返回值和诊断

### `V`

节点价值函数，大小为 `N`×`R`。

### `x`

最终控制策略。单 regime 时通常为 `N`×`p`；多个 regime 且控制列数一致时为 `N`×`p`×`R`；若各 regime 控制维数不同，则返回 cell 数组。

### `resid`

节点上的最终残差，大小为 `N`×`R`。普通 HJB 是 Hamiltonian 减贴现价值；有停止或干预时是完整 VI/QVI 残差。

### `info`

v4 返回标量结构体，常用字段包括：

```matlab
info.converged
info.iterations
info.residual_max
info.error_history
info.residual_history
info.grid
info.branches
info.adversary
info.stop_mask
info.intervention_mask
info.obstacle
info.target
info.destination
info.intervention_choice
info.intervention_action
info.lcp
```

其中：

- `stop_mask`：停止障碍活动的节点；
- `intervention_mask`：脉冲干预活动的节点；
- `branches`：每个节点选中的控制候选编号；
- `adversary`：每个节点选中的 Nature 控制（若有）；
- `lcp.residual`：LCP 后端返回的互补条件残差。

### `A`

最终政策下的稀疏生成矩阵，大小为 `N*R`×`N*R`。稳健模型中它对应求解器选中的最坏情形政策；它不是自动生成的物理分布方程矩阵。

## 8. 二维人力资本例子

```matlab
run_setup;
[model,amin,amax,hmin,hmax] = example_model_humancapital();

F = fundefn('lin',[36 23],[amin hmin],[amax hmax]);
nodes = funnode(F);
p = model.params{1};
G = hjb_grid(F,nodes,'fd');
a = G.S(:,1); h = G.S(:,2);

ell0 = min(1,max(0,((p.delta/p.theta).*h.^(1-p.alpha)).^(1/p.alpha)));
c0 = p.r*a + p.w*h.*(1-ell0);
V0 = c0.^(1-p.gamma)/(1-p.gamma)/p.rho;

opts = struct('maxit',300,'Delta',100,'tol',1e-6,...
    'res_tol',2e-5,'showiters',true);
[V,policy,resid,info] = HJBsolve(model,F,nodes,V0,opts);

fprintf('converged=%d, iterations=%d, residual=%.3e\n',...
    info.converged,info.iterations,info.residual_max);
figure;
subplot(1,2,1); imagesc(nodes{1},nodes{2},reshape(V,G.n)');
axis xy; colorbar; xlabel('Assets'); ylabel('Human capital');
subplot(1,2,2); imagesc(nodes{1},nodes{2},reshape(policy(:,2),G.n)');
axis xy; colorbar; xlabel('Assets'); ylabel('Human capital');
```

该示例使用四个上风方向组合，并在资产和人力资本边界上显式处理可行控制。

## 9. 常见错误

### `hjb:V0 Initial V must be N-by-regimes`

`V0` 的元素总数必须等于 `G.N * model.nregimes`。多 regime 时不要只传一个 `N`×1 向量。

### `hjb:grid Grid does not match fspace.n`

`snodes` 的每一维节点数必须与 `fspace.n` 一致，且节点严格递增。

### `hjb:boundary Outward drift at state-constraint boundary`

某个候选在状态约束边界产生了向域外漂移。应在候选回调中将该候选标记为不可行，或补充一个满足边界约束的候选。

### `hjb:generator Invalid Markov generator`

检查漂移、扩散、交叉扩散和边界类型。生成矩阵要求行和接近零，非对角元非负。

### `hjb:lcp` 或 `hjb_lcp_solve` 相关错误

确认 `lcpsolve` 在路径中；检查障碍值、初值和 `info.lcp.residual`。可以尝试 `opts.lcp_backend='moll'`，但需要另外安装兼容的 `LCP.m`。

### `hjb:collocation_boundary`

这是 `scsolve_jump` 的限制：配点法目前不支持“障碍问题 + 反射 Neumann 边界”的组合。此限制不适用于普通 FD 的 `HJBsolve`。

### 不收敛但数值看起来合理

先查看 `info.error_history` 和 `info.residual_history`，再进行网格加密。对于非光滑策略、自由边界或非凸问题，优先使用 FD 结果作为基准，不要只看节点残差。

## 10. 推荐工作流

```matlab
% 1. 确认路径
run_setup;
which HJBsolve -all;

% 2. 建立线性 FD 网格和初值
F = fundefn('lin',N,lower,upper);
nodes = funnode(F);

% 3. 运行 HJBsolve
[V,x,resid,info,A] = HJBsolve(model,F,nodes,V0,opts);

% 4. 检查收敛和生成矩阵
assert(info.converged);
assert(info.residual_max < opts.res_tol);
hjb_check_generator(A,opts.generator_tol,false);

% 5. 加密网格并比较 V、x 和残差
```

最后一步是必要的数值检查，尤其适用于非凸生产技术、停止边界和脉冲干预问题。
