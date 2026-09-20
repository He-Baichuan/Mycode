---
title: ABB(2017) 方法与代码实现梳理
tags: [结构估计, 收入动态, 消费, 分位数回归, 面板数据]
date: 2026-09-20
---

# Arellano–Blundell–Bonhomme (2017) 方法梳理与代码实现总结

> **论文**：Arellano, M., Blundell, R., & Bonhomme, S. (2017). *Earnings and Consumption Dynamics: A Nonlinear Panel Data Framework*. **Econometrica**, 85(3), 693–734. DOI: [10.3982/ECTA13795](https://doi.org/10.3982/ECTA13795)
>
> **代码**：`./code/Codes_ABB/`，25 个 `.m` 文件（约 17,800 行）+ 2 套 R 生命周期求解器。
>
> **本文档的依据**：论文正文已逐节阅读（§2 模型、§3 消费规则、§4 识别、§5 估计策略、§6 实证结果）；代码为逐文件阅读 + 直接解析 `.mat` 验证。凡涉及代码的论断都标注了 `文件:行号`；凡属本文档作者的独立复算，均单独标注。论文公式序号沿用原文编号，便于对照。

---

## 0. 一句话概括

ABB 用一个**分位数形式的非线性 Markov 过程**替代传统的线性永久/暂时收入分解，从而允许"冲击的持续性依赖于冲击的大小、符号和家庭收入史"；再用一个**非线性的消费规则**把收入的两个成分映射到消费，从而得到随资产、冲击水平和收入史而变的**部分保险度量**。估计上采用**随机 EM**：E 步用随机游走 MH 抽取潜在持续成分 $\eta$，M 步给定 $\eta$ 抽值做**分位数回归**。

---

## 1. 论文定位：相对 BPP(2008) 的突破

线性收入模型（随机游走永久/暂时模型、AR(1) 等）有一个共同性质：**所有冲击对应同一个持续性**，与家庭的收入史无关。ABB 的切入点是：这个性质是**假设**，不是事实。

论文的目标不是拒绝线性模型本身，而是**在不预设函数形式的前提下度量持续性**，并进而考察这一非线性如何传导到消费。为此他们做到：

| 维度 | 线性模型（BPP 2008 等） | ABB(2017) |
|---|---|---|
| 收入过程 | 永久 + 暂时 / AR(1) | 一般一阶 Markov，用条件分位数刻画 |
| 持续性 | 常数 | $\rho_t(\eta_{i,t-1},\tau)$，随状态**和冲击分位**变化 |
| 消费规则 | 消费**增长率**的线性投影 | 消费**水平**的非线性函数 $g_t(a,\eta,\varepsilon,\nu)$ |
| 保险度量 | 常数 $\phi,\psi$ | $\phi_t(a,\eta,\varepsilon)$、$\psi_t(a,\eta,\varepsilon)$ |
| 估计 | 矩条件 | 随机 EM + 分位数回归 |

> [!important] 一个关键的技术细节
> $\rho_t(\eta_{i,t-1},\tau)$ **可以取负值，绝对值也可以超过 1**（原文脚注 12）。这不是病态——若 $\ln\eta_{it}$ 是带高斯新息的标准随机游走，则 $\eta_{it}$ 本身是乘性随机游走，其持续性在中位冲击处为 1，在下四分位处约为 0.5，在上四分位处约为 2。所以"持续性"在这里不是自回归系数，而是**条件分位数函数的斜率**。

---

## 2. 模型

### 2.1 收入过程（原文 §2，式 (1)–(2)）

$$\boxed{y_{it} = \eta_{it} + \varepsilon_{it}}, \qquad i=1\dots N,\ t=1\dots T \tag{1}$$

- $y_{it}$：**税前劳动收入的对数**，已剔除一整套年龄虚拟变量；
- $\eta_{it}$：**持续成分**，服从一般一阶 Markov 过程；
- $\varepsilon_{it}$：**暂时成分**，零均值、跨期独立、与所有 $\eta_{is}$ 独立。

不去参数化 $\eta$ 的条件**分布**，而是用条件**分位数**表示：

$$\eta_{it} = Q_t(\eta_{i,t-1}, u_{it}), \qquad (u_{it}\mid \eta_{i,t-1},\eta_{i,t-2},\dots)\sim \text{Uniform}(0,1), \quad t=2\dots T \tag{2}$$

式 (2) 是**不失一般性的**——对任意连续条件分布都成立。$\eta$ 的依赖结构除了"一阶 Markov"之外没有任何限制。

> [!note] 原文明确交代的三点
> 1. **测量误差**：*"在缺乏额外信息的情况下，无法把暂时新息与经典测量误差区分开。因此我们估计的 $\varepsilon_{it}$ 分布，其解释是**暂时冲击与测量误差的混合**。"*（§2.1）——这一点常被引用时忽略。
> 2. **年龄**：$\eta$ 与 $\varepsilon$ 都假设**均值独立于年龄**，但条件分位数函数 $Q_t$ 与 $\varepsilon_{it}$ 的边际分布**都可以依赖 $t$**。由于同队列内年龄与日历时间完全共线，这种依赖同时吸收了年龄效应和总量冲击。
> 3. **初期条件**：$\eta_{i1}$ 的分布不作任何限制。

**canonical（经典）模型**是式 (3)：$\eta_{it} = \eta_{i,t-1}+v_{it}$，即 $\eta$ 为随机游走。此时 $\rho_t \equiv 1$，与状态和冲击无关。

### 2.2 非线性持续性、异方差与偏度（原文 §2.2，式 (4)–(6)）

$$\rho_t(\eta_{i,t-1},\tau) = \frac{\partial Q_t(\eta_{i,t-1},\tau)}{\partial \eta}, \qquad \rho_t(\tau) = \mathbb{E}\!\left[\frac{\partial Q_t(\eta_{i,t-1},\tau)}{\partial \eta}\right] \tag{4}$$

$\rho_t(\eta_{i,t-1},\tau)$ 度量"当 $\eta_{i,t-1}$ 被一个**秩为 $\tau$ 的当期冲击**击中时，它的持续性"。它同时依赖滞后成分和冲击分位——这正是"持续性依赖冲击大小与符号"的形式化。ABB 称之为**收入史的持续性**（persistence of earnings histories）。

**条件异方差**（原文 §2.2）：

$$\sigma_t(\eta_{i,t-1},\tau) = Q_t(\eta_{i,t-1},\tau) - Q_t(\eta_{i,t-1},1-\tau)$$

**条件偏度**（式 (6)，Kim–White 2004 型）：

$$\text{sk}_t(\eta_{i,t-1},\tau) = \frac{Q_t(\eta_{i,t-1},\tau)+Q_t(\eta_{i,t-1},1-\tau)-2Q_t(\eta_{i,t-1},\tfrac12)}{Q_t(\eta_{i,t-1},\tau)-Q_t(\eta_{i,t-1},1-\tau)} \tag{6}$$

**条件峰度**（脚注 14）：$\text{kurt} = \dfrac{Q(1-\alpha)-Q(\alpha)}{Q(\tau)-Q(1-\tau)}$。

> [!check] 这几条解决了代码里的一个"疑似笔误"
> `Code_earnings_fit.m:523` 的 `var2=MM(:,Ntau)-MM(:,1)` 被注释为 "conditional variance"，我一开始怀疑是命名错误。**实际上它完全正确**——这正是原文 §2.2 定义的 $\sigma_t$（分位差形式的离散度），不是方差。同样，`skew2` 与 `kurt1` 也严格对应式 (6) 和脚注 14。代码与论文在这里是精确对应的。

### 2.3 消费规则（原文 §3.2，式 (10)–(12)）

$$c_{it} = g_t(a_{it}, \eta_{it}, \varepsilon_{it}, \nu_{it}), \qquad t=1\dots T \tag{10}$$

- $a_{it}$：**对数资产**（同样剔除年龄效应）；
- $\nu_{it}$：消费函数中**未被观测的自变量**——对 $\nu$ 单调，跨期独立且独立于 $(a_{it},\eta_{it},\varepsilon_{it})$；
- $\nu$ 的经济解释是**偏好冲击**（提高边际效用）；**不失一般性把 $\nu_{it}$ 的边际分布标准化为标准均匀分布**；
- 原文明确指出：*"$\nu_{it}$ 的存在也可能部分吸收了消费支出中的测量误差。"*

**资产被设定为序贯外生（predetermined）**（式 (11)）：

$$a_{it} = h_t(a_{i,t-1}, c_{i,t-1}, y_{i,t-1}, \eta_{i,t-1}, \upsilon_{it}) \tag{11}$$

$h_t$ 对最后一个变量递增，$\upsilon_{it}$ 为 i.i.d. 均匀。标准线性资产规则（式 (7)，即预算约束）是 (11) 的**特例**，因此估计时**不需要**施加 (7)；但**模拟生命周期冲击时必须**对预算约束表态。

**部分保险度量**（式 (12)）：

$$\phi_t(a,\eta,\varepsilon) = \mathbb{E}\!\left[\frac{\partial g_t(a,\eta,\varepsilon,\nu_{it})}{\partial \eta}\right] \tag{12}$$

脚注 17 给出对暂时成分的对应量 $\psi_t(a,\eta,\varepsilon)=\mathbb{E}[\partial g_t/\partial\varepsilon]$。$1-\phi_t(a)$ 即"对持续冲击的保险程度"。

**冲击的动态效应**（这是脉冲响应的理论基础）：由链式法则，

$$\mathbb{E}\!\left[\frac{\partial}{\partial u}\Big|_{u=\tau} g_t\big(a, Q_t(\eta,u),\varepsilon,\nu_{it}\big)\right] = \phi_t\big(a,Q_t(\eta,\tau),\varepsilon\big)\cdot\frac{\partial Q_t(\eta,\tau)}{\partial u}$$

且 $\partial Q_t(\eta,\tau)/\partial u$（冲击对收入的即期反应）**对 $\eta$ 的导数等于 $\partial\rho_t(\eta,\tau)/\partial\tau$**——这把持续性度量和冲击反应直接联系起来了。原文：*"在实证分析中，我们将报告这些导数效应的**有限差分**对应物（'脉冲响应'）。"*

### 2.4 主要实证发现（原文 §6）

**收入侧**：

- 持续性最高的情形：**高收入家庭遇好冲击**（高 $\tau_{\text{init}}$、高 $\tau_{\text{shock}}$）、**低收入家庭遇坏冲击**，两者持续性都接近 **0.9–1**；
- 持续性最低的情形：**高收入遇坏冲击**、**低收入遇好冲击**，低至 **0.3–0.4**（Figure 1，PSID 与挪威数据一致）；
- 条件偏度：$\eta_{it}$ 在 $\eta_{i,t-1}$ 低时**右偏**、高时**左偏**；对数收入残差也有同向但不那么强的证据（Figure 4）；
- 对数收入增长**非高斯、负偏、高峰度**，存在 ARCH 效应（与 Guvenen et al. 2015 的行政数据一致）；
- **挪威数据**：2,873 户平衡子样本（2000–2005），非线性持续性与条件偏度pattern 与 PSID 一致；但**暂时成分 $\varepsilon$ 的离散度远小于 PSID**——暗示 PSID 中存在较大测量误差，或挪威的真实暂时冲击确实更小。

**消费侧**：

- 消费对收入的条件均值导数在 **0.2–0.3** 之间；**年长、高资产家庭的消费与收入变动关联更弱**；
- $\phi_t(a)$ 平均在 **0.3–0.4**，即*"超过一半的税前家庭收入波动被有效保险"*；年长且高资产家庭保险更充分。

**脉冲响应（§6.4，Figure 7–9）**：

- 收入：$\tau_{\text{shock}}=0.10$ 的负冲击使低收入家庭（$\tau_{\text{init}}=0.10$）收入下降 **7%**，而高收入家庭（$\tau_{\text{init}}=0.90$）下降 **19%**；
- 消费：同样一个负冲击使低收入家庭消费下降 **2%**，高收入家庭下降 **8%**；正冲击（$\tau_{\text{shock}}=0.90$）使高收入家庭消费上升 **5%**，低收入家庭上升 **11%**；
- 按年龄与资产分组（Figure 9）：高收入家庭的负冲击在**更晚年龄**影响更大——53 岁受冲击时收入下降 **40%**，而 37 岁受冲击时约 **20%**；**资产持有能削弱对负冲击的消费反应**（尤其对生命周期后期受冲击的家庭），但对正冲击的保险程度影响不大；
- 原文脚注 40 提醒：所谓"大冲击"是**相对**说法，因为它对应的是不同条件分布的分位数。

---

## 3. 识别与估计

### 3.1 识别（原文 §4）

模型是非线性状态空间模型。识别思路承自 Hu–Schennach (2008)、Hu–Shum (2012) 的非经典测量误差/潜变量文献，在**条件独立**约束下建立非参数识别。要点：

- **收入过程**（§4.1）：利用 $\varepsilon$ 跨期独立、$\eta$ 一阶 Markov 的结构，把观测 $y$ 的联合分布"反卷积"出 $\eta$ 的转移分布与 $\varepsilon$ 的边际分布；
- **消费规则**（§4.2）：把 BPP(2008) 的 IV 论证推广到非参数情形，关键是 $a_{it}$ 只需**序贯外生**；
- **家庭不可观测异质性**（§4.3）：加入时不变因子 $\zeta_i$ 后的识别。
- 补充材料给出了 $T=3$ 情形的严格证明（基于条件期望算子 $L_{y_2|y_1}$、$L_{\eta_2|y_1}$ 与完备性假设）。

### 3.2 估计策略：随机 EM，但**不是**似然方法（原文 §5.2）

估计目标可以紧凑地写成：

$$\hat\theta = \arg\min_\theta \mathbb{E}\!\left[\int R(y_i,\eta;\theta)\, f_i(\eta;\theta)\,d\eta\right]$$

其中 $R$ 是已知函数，$f_i(\cdot;\theta)=f(\cdot\mid y_i^T,\text{age}_i^T;\theta)$ 是给定收入数据后 $(\eta_{i1},\dots,\eta_{iT})$ 的**后验密度**。

> [!important] 原文的关键说明
> 该算法*"与随机 EM 算法（Celeux and Diebolt 1993）密切相关"*，但有一个重要区别：*"与 EM 不同，我们的**问题不是基于似然的**。我们利用分位数回归的计算便利性，在每一步 M-step 中**用一串分位数回归替代似然最大化**。"*
>
> 这一点在代码里看得很清楚：`Code_earnings_est.m` 的 M 步是 11 个独立的 `fminunc(@wqregk_*, ...)`（每个 $\tau$ 一个），而不是一个联合的似然优化。

算法（从 $\theta^{(0)}$ 出发，迭代 $s=0,1,2,\dots$ 直到 $\theta^{(s)}$ 收敛）：

```
1. 随机 E 步：从 f_i(·; θ⁽ˢ⁾) 抽取 η_i⁽ᵐ⁾ = (η_i1⁽ᵐ⁾,...,η_iT⁽ᵐ⁾), m = 1..M
2. M 步：    θ⁽ˢ⁺¹⁾ = argmin_θ  Σ_i Σ_m R(y_i, η_i⁽ᵐ⁾; θ)
```

**M 步的具体形式**——以 $a^Q_{k\ell}$ 为例：

$$\min_{(a^Q_{0\ell},\dots,a^Q_{K\ell})}\ \sum_{i=1}^N\sum_{t=2}^T\sum_{m=1}^M \rho_{\tau_\ell}\!\left(\eta^{(m)}_{it} - \sum_{k=0}^K a^Q_{k\ell}\,\varphi_k\big(\eta^{(m)}_{i,t-1},\text{age}_{it}\big)\right),\quad \ell=1\dots L$$

其中 $\rho_\tau(u)=u(\tau-\mathbb{1}\{u\le 0\})$ 就是检验函数（**与 `wqregk_pt_age.m` 中 `Obj=mean((...).*(tau-(...<0)))` 完全一致**）。这是一组标准分位数回归，目标函数凸。

**E 步**：*"由于似然函数有闭式表达，E 步是直接的。实践中我们使用随机游走 Metropolis–Hastings 抽样器，目标接受率约 **30%**。"*

### 3.3 论文明确交代的实现参数（与代码逐条对照）

| 项目 | 论文 §5.2 原文 | 代码 | 一致？ |
|---|---|---|---|
| 抽样数 $M$ | $M=1$ | `Mdraws=1` | ✅ |
| 收入迭代数 $S$ | 500 次迭代 | `Code_earnings_est.m:70` `maxiter=500` | ✅ |
| 每次迭代 MH 抽样 | 200 次 | `:73` `draws=200` | ✅ |
| 消费迭代数 | 200 次迭代、每次 200 抽样 | `Code_consumption_est.m:125,128` `maxiter=200; draws=200;` | ✅ |
| 平均区间 $\bar S$ | $\bar S=S/2$ | `Resqfinal=mean(Resqnew(:,:,(maxiter/2):maxiter))` | ✅ |
| 年龄效应 | *"先对年龄的四次多项式回归"* | `Code_earnings_est.m:42-51`，4 阶 Hermite | ✅ |
| $\varepsilon$ 与年龄不相关 | *"每次迭代都施加"* | 归一化在 `:396-397` 每次迭代执行 | ✅ |
| 分位点网格 $L$ | $L=11$，$\tau_\ell=\ell/(L+1)$ | `Ntau=11`，`Vectau=(1/12:1/12:11/12)` | ✅ |
| 分位数函数形式 | $[\!\tau_1,\tau_L]$ 上**分段线性**；$a^Q_0$ 在两端为**指数分布分位数** | 代码的线性插值 + Laplace 尾 | ✅ |
| 基函数 | *"低阶 Hermite 多项式乘积，每项取标准化变量"* | `hermite.m` + `(x-mean)/std` | ✅ |
| 消费规则阶数 | 脚注 35：*"消费规则估计使用次数为 $(2,2,1)$ 的 Hermite 张量积"* | `M1=2,M2=2,M3=1,M4=1` | ✅ |

**关于起始值**（这条解释了一个容易误会的代码特征）：

> *"我们的实验中发现算法可能'卡'在马尔可夫链的某个局部状态。我们**从大量不同的初始参数值出发**运行算法，并选择在迭代中平均对数似然最高的那一组估计。"*

这正好解释了代码里的两件事：
1. `rng('shuffle')` 不是疏忽，而是**多起点探索**机制的一部分；
2. 每次迭代监控并打印 `mat_lik`（`Code_earnings_est.m:466`、`Code_consumption_est.m:542`），末尾打印 `mean(mat_lik((maxiter/2):maxiter))`（`:517` / 消费脚本同）——**这就是论文所说的选择准则**。

**序贯 vs 联合**（§5.2 开头）：

> *"我们选择序贯估计而非联合估计 $(\theta,\mu)$，理由是 $\theta$ 仅从收入过程即可识别。相反，在联合估计中，收入过程的估计会被消费模型部分驱动。"*

**渐近性质**（§5.2 "Properties"）：Nielsen (2000) 给出随机 EM 链的遍历性与渐近分布；Arellano–Bonhomme (2016) 在**非似然、基于分位数估计方程**的情形下刻画了 $\hat\theta$ 的渐近分布。在参数模型正确设定、$K$ 与 $L$ **固定**时，$\hat\theta$ 是 $\sqrt N$ 相合且渐近正态的。

### 3.4 归一化与指数尾的精确含义

原文脚注 24 给出 $a^Q_0$ 的完整表达式：

$$a^Q_0(\tau) = \frac{1}{\lambda^Q_-}\log\!\Big(\frac{\tau}{\tau_1}\Big)\mathbb{1}\{0<\tau<\tau_1\} + \sum_{\ell=1}^{L-1}\Big(a^Q_{k\ell} + \frac{a^Q_{k,\ell+1}-a^Q_{k\ell}}{\tau_{\ell+1}-\tau_\ell}(\tau-\tau_\ell)\Big)\mathbb{1}\{\tau_\ell\le\tau<\tau_{\ell+1}\} - \frac{1}{\lambda^Q_+}\log\!\Big(\frac{1-\tau}{1-\tau_L}\Big)\mathbb{1}\{\tau_L\le\tau<1\}$$

这与代码里的 Laplace 尾表达式逐项对应（`b1` $=\lambda^Q_-$，`bL` $=\lambda^Q_+$）：

```matlab
Mateta_true(:,1)=Mateta_true(:,1)+((1/(b1true_e0)*log(V_draw/Vectau(1))).*(V_draw<=Vectau(1))...
    -(1/bLtrue_e0*log((1-V_draw)/(1-Vectau(Ntau)))).*(V_draw>Vectau(Ntau)));
```

**归一化的含义**（这段是我基于上述公式推导的，代码只写了结果）：

`Code_earnings_est.m:396-397` 两条限制

```matlab
Resqnew_eps(:,:,iter)=Resqnew_eps(:,:,iter)-mean(Resqnew_eps(:,:,iter)')'*ones(1,Ntau);
Resqnew_eps(1,:,iter)=Resqnew_eps(1,:,iter)-((1-Vectau(Ntau))/bL_eps-Vectau(1)/b1_eps)*ones(1,Ntau);
```

- 第 1 条：把 Hermite 系数**减去系数均值**。由于 $\varphi_0\equiv 1$，这等价于把分布的**位置**归一化为零；
- 第 2 条：对第一个系数（截距）再平移 $-\big[(1-\tau_L)/\lambda_+ - \tau_1/\lambda_-\big]$。对指数尾求积分可知，左尾在 $(0,\tau_1]$ 上对均值的贡献是 $-\tau_1/\lambda_-$，右尾在 $[\tau_L,1)$ 上的贡献是 $(1-\tau_L)/\lambda_+$。因此这一平移恰好**扣除了两段指数尾对均值的净贡献**——两条合起来，就是把整个分布的均值精确地钉在零。

> [!warning] 一个如实记录的实现事实
> 基线 `Code_earnings_est.m` 中，这两条归一化**只施加于 $\varepsilon$**；$\eta$ 的转移系数 `Resqnew` 与初期条件 `Resqnew_e0` **没有**对应的归一化行（全文只有 :396-397 两行）。FE 版 `Code_earnings_est_FE.m:436-440` 同样只对 `eps` 与 `zeta` 归一化。这不是笔误的可能性很高——因为原文对 $\eta$ 的位置没有施加归一化，$\eta$ 的水平由初期条件的年龄基与 Laplace 参数共同吸收。但我无法从论文中确证这一点，故如实记录而不作断言。

---

## 4. 代码地图

### 4.1 主流程

| 文件 | 作用 | 输入 | 输出 | 论文图 |
|---|---|---|---|---|
| `first_stage_regs.do` / `select.do` | Stata 建样本、一阶段回归、导出残差 | `data.dta`, `select.dta` | `first_stage_regs.out` | — |
| `Code_earnings_est.m` | **收入模型估计** | `first_stage_regs.out` | `data_hermite2.mat` | — |
| `Code_earnings_fit.m` | 收入模型模拟与拟合检验 | `data_hermite2.mat` | — | 1, **2abd**, 3, 4，S1–S8 |
| `Code_consumption_est.m` | **消费模型估计** | `data_hermite2.mat` | `data_hermite_cons2.mat` | — |
| `Code_consumption_fit.m` | 消费模型模拟 | `data_hermite_cons2.mat` | — | 5, 6，S20ab, S21a |
| `Code_IR_cons.m` | 脉冲响应（两种资产规则） | `data_hermite_cons2.mat` | — | 7, 8，S31 |
| `Code_IR_cons_by_assets.m` | 按年龄与初始资产分组的脉冲响应 | 同上 | — | 9，S34(上) |
| `Code_canonical_model_est_and_fit.m` | 经典线性模型对照 | `data_hermite_cons2.mat` | — | **2c** |
| `Code_canonical_model_IR.m` | 线性模型的脉冲响应 | `data_hermite_cons2.mat` | — | 7, 8 的 (g)(h) |

### 4.2 辅助函数

| 文件 | 作用 |
|---|---|
| `hermite.m` | **概率论者 Hermite 多项式** $He_n(x)$，满足 $He_n' = n\,He_{n-1}$ |
| `rq.m` | Koenker 的 Frisch–Newton 内点法分位数回归（由 Ox 转写；`beta=0.9995, small=1e-5, max_it=50`） |
| `fun_qrlocal.m` | 局部线性分位数回归（**未被主流程调用**，遗留） |
| `postr_QRMCMC_age_hermite.m` | 收入模型的个体后验密度 $f_i(\eta;\theta)$ |
| `postr_nonlinear_consumption_age_hermite_predet.m` | 消费模型的个体后验密度 |
| `wqregk_e0_age.m` / `wqregk_pt_age.m` / `wqregk_eps_age.m` | 收入 M 步三个检验函数目标 |
| `covariance_cm.m` / `covariance_cm_consumption_twostep.m` | 经典模型的等权最小距离目标 |

### 4.3 扩展模块

| 文件 | 作用 | 论文图 |
|---|---|---|
| `Code_earnings_est_FE.m` / `_fit_FE.m` | 收入模型 + 不可观测异质性 $\zeta_i$ | S15 |
| `Code_consumption_est_FE.m` / `_fit_FE.m` | 消费模型 + 家庭固定效应 $\xi_i$ | S24, S25, S20c, S21b |
| `wqregk_zeta_age.m`, `wqregk_e0_age_FE.m`, `wqregk_eps_age_FE.m`, `postr_QRMCMC_age_hermite_FE.m`, `postr_nonlinear_consumption_age_hermite_FE_predet.m` | FE 版辅助函数 | — |
| `Code_earnings_{nonparametric,parametric}_bootstrap{,_fit}.m` | 收入的两种 bootstrap | S9–S14 |
| `Code_consumption_{nonparametric,parametric}_bootstrap{,_fit}.m` | 消费的两种 bootstrap | S22, S23 |
| `Code_IR_cons_nonparametric_bootstrap{,_fit}.m` | 脉冲响应的 bootstrap 置信带 | S26–S29 |
| `Code_IR_cons_FE.m` / `_by_assets.m` | FE 版脉冲响应 | S30, S32, S33, S34(下) |
| `Code_simulated_model.m` | 外部生命周期模型（年度校准版） | S35, S36 |
| `Code_simu_nonlinear.m` / `_cons.m` | 外部生命周期模型（双年估计版） | S37 |
| `Simulations/Calibrated/*.r`、`Simulations/Estimated/*.r` | R 生命周期求解器 + `run.cons.reg.st.parallel.r` | — |

---

## 5. 代码级拆解

### 5.1 数据构建（Stata → MATLAB）

`first_stage_regs.do`：

1. **平衡子样本**：`keep if numwav == 6`（每个体保留 6 期）；
2. **剔除缺失**：资产、消费、劳动收入、教育、州、出生年、孩子数、家庭规模；
3. **一阶段回归**：对 `log_cons`、`log_totly`、`log_ass` 分别
   ```
   xi i.educ*i.yb i.weduc*i.wyb i.state i.fsize i.kids i.race i.wrace
   reg <depvar> kidsout bigcity extra _I*
   predict u<depvar>, res
   ```
4. **异方差标准化**：残差平方对 `i.educ i.yb i.weduc i.wyb` 回归，取 $\hat u/\sqrt{\widehat{u^2}}$；
5. **导出**：`outsheet utoty uc ua age using first_stage_regs, nonames replace`。

**实测**：`first_stage_regs.out` 共 **4,752 行 = 792 个体 × 6 期**，列为 `utoty uc ua age`。这与原文 §6.1 的 PSID 样本（工作家庭，1999–2009 双年）一致。

MATLAB 侧（`Code_earnings_est.m:15-35`）读入展平向量，`T=6`、`N=size(Y,1)/T`，再 reshape 成 $N\times T$。**年龄效应**在 `:42-51` 用 4 阶 Hermite 年龄回归剔除（对应原文的"age quartic"）。消费与资产在消费脚本中同样处理（`:34`、`:52`）。

### 5.2 收入估计 `Code_earnings_est.m`

**参数设置**

| 参数 | 值 | 对应论文 |
|---|---|---|
| `N`, `T` | 792, 6 | §6.1 |
| `Ntau`, `Vectau` | 11, (1/12,…,11/12)′ | $L=11$，$\tau_\ell=\ell/(L+1)$ |
| `K1, K2` | 3, 2 | 式 (21) 的 Hermite 乘积（$4\times3=12$ 项） |
| `K3` | 2 | $\eta_{i1}$ 对年龄 |
| `K4` | 2 | $\varepsilon_{it}$ 对年龄 |
| `maxiter`, `draws`, `Mdraws` | 500, 200, 1 | $S=500$，200 MH draws，$M=1$ |
| `var_prop1..6` | .08,.03,.03,.03,.03,.05 | 逐期提议方差（论文未报告） |
| `b1, bL` | 10, 10 | $\lambda^Q_-,\lambda^Q_+$ 初值 |

**E 步的实现**：`postr_QRMCMC_age_hermite.m` 返回每个个体的完整数据后验密度（N×1，**非对数**），由两部分相乘：

| 成分 | 含义 | 分位数函数 |
|---|---|---|
| `denstot` | 数据似然：$\varepsilon=y-\eta$ 的密度 | `Resqinit_eps` |
| `dens2tot` | $\eta$ 的先验：初期条件 $\eta_1$ × 转移 $\eta_t\mid\eta_{t-1}$ | `Resqinit_e0`、`Resqinit` |

密度由分位数函数的分段线性插值隐含给出（区间 $j$ 上高度 $=(\tau_{j+1}-\tau_j)/(Q_{j+1}-Q_j)$），两端接指数尾。

**实现上的关键技巧——跨个体向量化**（`Code_earnings_est.m:255-330`）：MH 循环**一次同时更新全部 N 个个体**，每人独立链、独立接受/拒绝：

```matlab
Matdraw(:,1)=Nu_chain1(:,j-1)+sqrt(var_prop1)*randn(N,1);
newObj=postr_QRMCMC_age_hermite(Matdraw);
r=(min([ones(N,1) newObj./Obj_chain(:,j-1)]'))';
prob=rand(N,1);
Nu_chain1(:,j)=(prob<=r).*Matdraw(:,1)+(prob>r).*Nu_chain1(:,j-1);
```

$a \to b \to c$ 依次更新，即 component-wise MH within Gibbs。

### 5.3 消费估计 `Code_consumption_est.m`

> [!important] 这里我最初理解错了，读论文后纠正
> 我第一版的总结写"消费规则是 OLS 投影，不是分位数回归"，并把它当作一个建模选择来描述。**更准确的理解是**：消费规则**也是**分位数函数，只是被设定为**对 $\tau$ 可加**。原文式 (22)：
> $$g_t(a_{it},\eta_{it},\varepsilon_{it},\tau) = \sum_{k=1}^{K} b^g_k\,\varphi_k(a_{it},\eta_{it},\varepsilon_{it},\text{age}_{it}) + b^g_0(\tau)$$
> 并且由于*"数据没有显示消费偏离对数正态的证据"*，取 $b^g_0(\tau)=\alpha+\sigma\Phi^{-1}(\tau)$。
>
> **也就是说这是一个正态位置—尺度模型**：$\varphi_k$ 的系数 $(b^g_k)$ 与 $\tau$ **无关**，$\tau$ 只通过截距 $b^g_0$ 进入。于是在该设定下，条件分位数的系数向量与条件均值的系数向量**相同**，OLS 就是（有效）估计量，残差方差给出 $\sigma$。常数 $\alpha$ 被吸收进全零次数的 Hermite 项（$\varphi_0\equiv 1$）。
>
> 这与代码完全对应（`:470-484`）：

```matlab
% 消费规则基：(a, η, ε, age) 的 4 维张量积 Hermite，M=(2,2,1,1)
XX=[];
for mm1=0:M1, for mm2=0:M2, for mm3=0:M3, for mm4=0:M4
    XX=[XX hermite(mm1,(Atot_t-meanA)/stdA)...
           .*hermite(mm2,(Matdraw_t-meanY)/stdY)...            % η 的 MH 抽值
           .*hermite(mm3,(Ytot_t-Matdraw_t-meanY)/stdY)...     % ε = y − η
           .*hermite(mm4,(AGEtot_t-meanAGE)/stdAGE)];
end, end, end, end
Resnew(1:36,iter) = pinv(XX)*Ctot_t;                          % 即 b^g_k
Resnew(37,iter)   = mean((Ctot_t-XX*Resnew(1:36,iter)).^2);   % 即 σ²
```

**资产规则**（式 (24)，`:489-506`）：$a_{it}$ 对 $(a_{i,t-1},c_{i,t-1},y_{i,t-1},\eta_{i,t-1},\text{age}_{it})$ 的 Hermite 张量积回归，$R=(1,1,1,1,1)$ → $2^5=32$ 系数 + 1 方差，同样是"对 $\tau$ 可加"+ OLS。原文脚注 23 提到：早期版本曾施加"$\eta_{i,t-1}$ 不进入式 (24)"，结果定性相似。

**期初资产**（式 (23)，`:510-523`）：$a_{i1}$ 对 $(\eta_{i1},\text{age}_{i1})$ 的**分位数回归**，$M_5=M_6=1$ → 4 系数 × 11 分位点 + 指数尾。

**消费后验**（`postr_nonlinear_consumption_age_hermite_predet.m`）由 **5 项**相乘：

| 项 | 内容 | 形式 | 对应 |
|---|---|---|---|
| `denstot` | 收入数据似然（$\varepsilon=y-\eta$） | 分位数隐含密度 | 式 (2) |
| `dens2tot` | $\eta$ 先验（初值 × 转移） | 分位数隐含密度 | 式 (21) |
| `dens3tot` | 消费数据似然 | **正态**，方差 `Resinit(37)` | 式 (22) |
| `dens4tot` | 期初资产似然 | 分位数隐含密度 | 式 (23) |
| `dens5tot` | 资产规则似然 | **正态**，方差 `Resinit_a2(end)` | 式 (24) |

> [!note] 关于测量误差
> **收入侧**：原文明确说明 $\varepsilon_{it}$ 是"暂时冲击与测量误差的混合"，二者无法分离（§2.1）。
> **消费侧**：原文同样明确说明偏好冲击 $\nu_{it}$ *"也可能部分吸收了消费支出中的测量误差"*（§3.2）。
> 所以在代码里消费的随机成分全部由 `dens3tot` 的单一正态方差吸收，是**忠实于论文设定的**，不是代码简化。

### 5.4 E 步的初值处理

消费脚本在 MH 链启动前没有 $\eta$ 抽值，因此初值用观测收入的**随机对称拆分**作为 $(\eta,\varepsilon)$ 的粗糙代理（`Code_consumption_est.m:216-218`）：

```matlab
U=randn(N*T,1);
Y_t1=Y_t+U/2;    Y_t2=Y_t-U/2;    % 作为 (η, ε) 的代理进入初始 OLS
```

第一次迭代之后，$\eta$ 换为 MH 抽值、$\varepsilon$ 换为 $y-\hat\eta$。论文未描述这一初值细节（属于实现层面）。

### 5.5 拟合与图形 `Code_*_fit.m`

这些脚本**不估计任何参数**，是给定估计值的**参数化模拟 + 拟合检验**。流程（以 `Code_earnings_fit.m` 为例）：

1. 载入 `data_hermite2.mat`，估计量改名为 `*true`；
2. **扩样本**：`Nsim=20` → $N=792\times20=15840$；
3. **逆 CDF 抽样**：抽 $U\sim U(0,1)$，在 11 个估计分位点间**分段线性插值**反演，两端接解析指数尾（代码见 §3.4）；
4. **合成**：`Ytilde = Mateta_true + Mateps_true`；
5. **统计量**：用 `rq` + Hermite 导数法计算持续性曲面
   $$\texttt{Mat3*ResP\_data}=\frac{\partial Q_\tau(y_t\mid y_{t-1})}{\partial y_{t-1}}$$
   横轴 = $\tau_{\text{shock}}$，纵轴 = $\tau_{\text{init}}$；
6. **偏度/离散度/峰度**用 §2.2 的 Bowley 型分位统计量。

> [!note] 三张持续性图的基函数口径（初读容易误会）
> | 图 | 对象 | 基函数 | 列数 |
> |---|---|---|---|
> | Fig 2a（= Fig 1a/1b/S1 同源） | **数据**中 $y_t\mid y_{t-1}$ 的持续性 | `Mat1`：仅 $y_{t-1}$ 的 Hermite | $K_1{+}1=4$ |
> | Fig 2b | **模拟数据**中 $y_t\mid y_{t-1}$ 的持续性 | `MatS1`：同上 | $K_1{+}1=4$ |
> | Fig 2d / S3 | **模型隐含的 $\eta$** 的持续性 | `Resqtrue`：$\eta_{t-1}\times\text{age}$ | $(K_1{+}1)(K_2{+}1)=12$ |
>
> Fig 2a 与 2b 用**同一个 4 列基**，因此数据与模拟**可比**（这正是拟合检验的意义）；Fig 2d 画的是另一个对象（潜在 $\eta$ 而非观测 $y$），用含年龄交互的 12 列基是恰当的。原文 Figure 2 的注释也明确区分了这三者。

### 5.6 经典（canonical）线性对照模型

原文脚注 42 明确说明：canonical 模型的消费被设定为 $c_{it}$ 是 $\eta_{it}$、$\varepsilon_{it}$ 与一个跨期独立可加误差的线性函数，*"该模型用基于协方差约束的**等权最小距离**估计"*。这正好就是 `covariance_cm.m` / `covariance_cm_consumption_twostep.m` 的做法。

**第一步**：6×6 样本协方差对 3 个参数做等权最小距离

$$\text{Sig}(t,t)=\sigma^2_{\eta_1}+(t-1)\sigma^2_v+\sigma^2_\varepsilon,\qquad \text{Sig}(t,t')=\sigma^2_{\eta_1}+(\min(t,t')-1)\sigma^2_v$$

**第二步**：把 $(Y,C)$ 拼成 12×12 协方差矩阵，拟合 $\phi_\eta,\phi_\varepsilon,\sigma^2_\xi$（代码见 `Code_canonical_model_IR.m:67-119`）。

> [!check] 本文档作者的独立复算
> 用与上述目标**完全相同**的等权最小距离、3,000 次多起点全局优化：
>
> **第一步**：参数（以标准差计）$=(0.33598,\ 0.13478,\ 0.29896)$，目标值 $=0.009403$
> → $\sigma^2_{\eta_1}=0.11288$、$\sigma^2_v=0.01816$、$\sigma^2_\varepsilon=0.08938$（$Y$ 对角元 0.2106–0.3197，拟合良好）
>
> **第二步**：$\phi_\eta=0.6106$、$\phi_\varepsilon=\mathbf{-0.1793}$、$\sigma^2_\xi=0.06305$，目标值 $=0.038055$
>
> **关于 $\phi_\varepsilon$ 为负**：负号是**目标函数的性质，不是优化器陷入局部解**（3,000 次重启一致收敛到同一点）。合理的解读是：用 3 个参数去拟合完整的 12×12 协方差矩阵（144 个矩），线性模型无法匹配 (Y,C) 的协方差结构。这本身是 ABB 论点的一个侧证。原文正文并未报告 $\phi_\varepsilon$ 的具体数值，因此**建议按"论文未报告"处理，不要当作论文的正式估计引用**。

### 5.7 脉冲响应

**论文的定义**（§6.4）：在**年龄 35** 时所有家庭处于同一个持续性分位 $\tau_{\text{init}}$，在**年龄 37** 受到冲击 $\tau_{\text{shock}}$，然后报告**两条路径的年龄别中位数之差**（相对受 $\tau=0.50$ 中位冲击的家庭）。*"我们报告 100,000 次模型模拟的年龄别中位数。"* 原文脚注 41 说明：*"我们沿用'脉冲响应'这一说法"*，并引用 Gallant–Rossi–Tauchen (1993)、Koop–Pesaran–Potter (1996) 关于非线性模型中 IRF 的工作。

**代码实现**（`Code_IR_cons.m`）与论文逐项对应：

```matlab
N=100000;            % 原文：100,000 simulations
aa_ref=35; nage=13;  % 起始年龄 35，覆盖 35–59
tau_init=.1;  tau_shock=.9;

% 第一次模拟（基准）
tau0=.50;
...
if jj==1, V_draw=tau0*ones(N,1); else V_draw=unifrnd(0,1,N,1); end   % 仅年龄 37 那期固定
% 第二次模拟（处理）
tau0=tau_shock;
...
ResIR = nanmedian(Y) - nanmedian(Y1);    % 收入响应
ResIR = nanmedian(C) - nanmedian(C1);    % 消费响应
```

即两条路径**除了 37 岁那一期的冲击分位之外完全相同**（共享 `rng` 状态），因此逐个体可比。

**两种资产规则**（`acc` 开关，`Code_IR_cons.m:224-252`）：

| `acc` | 规则 | 对应 | 代码 |
|---|---|---|---|
| `2`（论文 Figure 7/8 所用） | 线性预算约束 + 3% 双年利率 + 资产下限 | 原文式 (7) | `A(:,jj+1)=log(max((1+.03)*exp(A(:,jj))+exp(Y(:,jj))-exp(C(:,jj)),floor_par*ones(N,1)));` |
| `1` | **从数据估计**的资产规则（32 系数 + 正态残差） | 原文式 (11) | `A(:,jj+1)=XXA*Restrue_a2(...)+sqrt(...)*randn(N,1);` |

> [!warning] 一个容易搞反的命名
> 代码注释把 `acc=1` 叫做 "nonlinear assets rule"，把 `acc=2` 叫做 "standard linear assets accumulation rule"。但从论文的角度看：
> - `acc=2` 是**理论上的预算约束**式 (7)（机械的会计恒等式，本身是线性的）；
> - `acc=1` 是**数据估计出来的非线性规则**式 (11)。
>
> 二者不是"对/错"的关系，而是论文刻意对比的两种设定。原文脚注 43：用估计的非线性资产规则（式 (11)）得到的结果*"与基准设定相比没有显著差别"*，对应 Figure S31。脚注 44 指出 Figure S34 的结果*"与 Figure 9 相比有一些差异，尤其是对正向收入冲击的响应"*。

**分组版本** `Code_IR_cons_by_assets.m`：结构相同，通过 `aa_ref`/`nage` 切换年轻组（35）与年老组（51），并按初始资产分位（0.10 虚线 / 0.90 实线）分组——正是原文 Figure 9 的注：*"初始资产在年龄 35（年轻家庭）或 51（年老家庭），位于 0.10 分位（虚线）和 0.90 分位（实线）。"*

**线性模型对照** `Code_canonical_model_IR.m`：$\eta$ 走随机游走，$c_t=\phi_\eta\eta_t+\phi_\varepsilon\varepsilon_t+\sqrt{\sigma^2_\xi}\xi_t$，用**完全相同**的有限差分方法算 IR，从而与非线性模型可比（对应原文 Figure 7/8 的 (g)(h) 面板）。

### 5.8 不可观测异质性（FE 扩展，原文 §5.1 末 + 式 (25)–(27)）

**消费侧**：原文式 (25)–(27)

$$c_{it} = g(a_{it},\eta_{it},\varepsilon_{it},\text{age}_{it},\xi_i,\nu_{it}) \tag{25}$$

**尺度归一化**（式 (26)）：

$$\sum_{k=1}^{K} b^g_k\,\varphi_k(0,0,0,\overline{\text{age}},\xi) = \xi \quad\text{对所有 }\xi$$

$$\xi_i = q(a_{i1},\eta_{i1},\text{age}_{i1},\omega_i),\qquad \omega_i\sim U(0,1) \tag{见 §5.1 末}$$

这正是代码 `Code_consumption_est_FE.m:548-571` 用**约束最小二乘**实现的：`r(2)=1` 把 $\xi$ 载荷的第一个系数**固定为 1**，`Vect_norm` 是各 Hermite 基在 0 点的取值乘积——与式 (26) 完全一致。消费后验相应地由 5 项变为 **6 项**（新增 `dens6tot`）。

注意：代码里 $\xi$ 的分位数函数用的是 $(A_1, Y_1, \text{AGE}_1)$ 而论文写的是 $(a_{i1},\eta_{i1},\text{age}_{i1})$——因为初期 $\eta_1$ 未观测，代码用观测 $Y_1$ 作代理（与 §5.4 的初值处理同源）。

**收入侧**：原文脚注 26——在 $y_{it}=\eta_{it}+\zeta_i+\varepsilon_{it}$ 中，*"我们通过另一个序列分位数模型，允许 $\eta_{i1}$、$\zeta_i$ 与 $\text{age}_{i1}$ 之间存在灵活的依赖"*。这正对应 `wqregk_e0_age_FE.m` 使用的 $\zeta\times\text{age}_1$ 张量积基。$\zeta$ 作为 MH 状态的**第 $T+1$ 列**（`Code_earnings_est_FE.m:337`），只做水平位移、不进入 $\eta$ 的动态。

**论文对 FE 结果的说法**：*"与 Figure 2 相比，允许家庭效应**降低了持续性**。而且非线性pattern 比同质情形**更明显**。当 $\tau_{\text{init}}$ 与 $\tau_{\text{shock}}$ 接近时持续性接近 1，但当大的正（负）冲击击中低（高）收入家庭时显著更低。"*（§6.2）消费侧：*"有异质性的模型中消费对收入变动的响应更小……且消费效应似乎更快回归中位数。"*（§6.4）

### 5.9 Bootstrap

**论文的定位**（§6.2 "Confidence Bands"）：比较两种方法。

> *"第一种是**在家庭层面聚类的非参数 bootstrap**。第二种是**参数 bootstrap**。后者要求参数模型正确设定，前者在误设定下仍可能相合且允许无限制的序列相关（**然而我们不知道在这个设定下它有正式的理论依据**）。"*

原文还说明：*"图 S9–S14 报告两种方法下的点态 95% 置信带，同时也给出基于非参数 bootstrap 的**一致（uniform）置信带**。所有情形下关于非线性持续性和条件偏度的主要发现都相当精确。同时，一致置信带更宽，对 $\eta$ 成分尤其如此。"*

**代码实现**：三种 bootstrap 脚本共享同一骨架：**(i)** 把主估计当作"真值"；**(ii)** 每个 replication **完整重估一遍模型**；**(iii)** 在重抽/模拟数据上算 fit 统计量并存盘。

| | 非参数 | 参数 |
|---|---|---|
| 重抽对象 | **个体**（整条时间序列，有放回）：`RDraw=floor(N*rand(N,1))+1; Y=Y(RDraw,:);` | **不重抽数据**，从估计模型生成全新面板 |
| 年龄曲线 | 每次重抽后**重新估计** | 固定 |
| replication 数 | 收入 500 / 消费 200 / IR 200 | 收入 500 / 消费 20 |

**共同点（重要）**：**都不从渐近分布抽参数**。参数不确定性完全通过"每次 replication 重估整个模型"来体现。

**置信带算法**（所有 `*_fit.m` 共用）：

1. 用某统计量的第一个 $\tau$ 判断 replication 是否成功：`ind=(Vect~=0).*(isnan(Vect)==0)`；
2. **点态 95% 带**：`quantile(Vect,.025)` / `quantile(Vect,.975)`；IR 的带做了 **recentered** 处理；
3. **一致带**：以点态带中点为心、按因子 $\lambda$ 放大
   ```matlab
   lambda=1.76;
   Res_fit_down2=(down+up)/2-lambda*(up-down)/2;
   ```
   $\lambda$ 取值：收入非参数 1.76/1.60/1.69/1.26/1.32；消费非参数 1.45/1.57；
4. **覆盖率校准**：`sum(prod.*(Vect>=down2).*(Vect<=up2))/Nboot`，即整条曲线都落在带内的比例。

**退化 replication**：接受率 > 0.95（链几乎不动）或似然为 NaN 的 replication 被丢弃并重抽。

### 5.10 外部生命周期模型（R）

`Simulations/` 下的 R 代码求解一个**标准生命周期消费–储蓄模型**，用来在一般生命周期框架中模拟 ABB 的估计结果。两套：

- **Calibrated**（年度，25–94 岁，`Twork=35, Tret=35`）：收入用"简单转换过程"；
- **Estimated**（**双年**，25–93 岁，`Twork=18, Tret=17`）：收入用**数据估计出的 ABB 非线性过程**。

模型要素：Gouveia–Strauss 税函数、养老金（`pencapfrac=2.2`）、自然借贷上限的倒推递归、完美年金市场（Estimated 版）、CRRA $\gamma=2$、$R=1.03$（年度）/ $1.03^2$（双年）、$\beta=1/R$。求解用倒推法 + `mclapply` 在 $z$ 网格上并行。

**文件命名**：

| 记号 | 含义 |
|---|---|
| `nl` | **n**on**l**inear：ABB 估计出的非线性收入过程 |
| `rw` | **r**andom **w**alk：经典线性过程 |
| `nbl` / `zbl` | **n**atural / **z**ero **b**orrowing **l**imit |
| `eps80` | 暂时冲击网格 `ngpe=80` |
| `samey` | 重新缩放使平均收入与 nl 一致 |

**产出链**：`run.model.*.r`（前向模拟，存 `.dat`）→ `R_matlab.r`（`load` + `writeMat`）→ `.mat` → MATLAB 绘图。`zbl.mat` / `nbl.mat` 各只含一个 **11×11 的 `persis` 矩阵**——消费响应对（资产分位 × 年龄分位）的曲面，由 `Simulations/Estimated/run.cons.reg.st.parallel.r` 制造。

---

## 6. 复现指南

### 6.1 运行顺序

```
Stata:  select.do ──▶ data.dta
        first_stage_regs.do ──▶ first_stage_regs.out
                    │
MATLAB:  Code_earnings_est.m ──▶ data_hermite2.mat
                    │
        ┌───────────┴────────────┐
        ▼                        ▼
 Code_earnings_fit.m      Code_consumption_est.m ──▶ data_hermite_cons2.mat
 (Figs 1,2abd,3,4,S1-S8)         │
                    ┌────────────┼────────────┬──────────────┐
                    ▼            ▼            ▼              ▼
        Code_consumption_fit  Code_IR_cons  Code_IR_cons_   Code_canonical_
        (Figs 5,6,S20ab,S21a) (Figs 7,8)    by_assets(Fig 9) model_*(Fig 2c)
```

### 6.2 已知的坑与注意事项

> [!warning] 运行前必读
> 1. **所有 `save` 语句在仓库中都被注释掉**——如 `Code_earnings_est.m:519` 的 `%save data_hermite2.mat`、`Code_consumption_est.m:588`、`Code_earnings_est_FE.m:578` 等。仓库直接提供预计算好的 `.mat`。若要重跑估计，需手动取消注释。
> 2. **bootstrap 脚本从硬编码的 Windows 路径读数据**——如 `Code_earnings_nonparametric_bootstrap.m:61` 的 `C:\Dfile\Dossiers\matlab\Qpanel\...\first_stage_regs.out`。在 macOS/Linux 上必须改写。
> 3. **`rng('shuffle')` 使结果不可复现**——该调用出现在 `Code_earnings_est.m:67`、`Code_consumption_est.m:122`、`Code_earnings_fit.m:69` 以及所有 bootstrap 和 IR 脚本中。**但请注意**：根据原文 §5.2，作者本来就是**从大量不同初值出发、按平均对数似然择优**（见 §3.3）。所以 `rng('shuffle')` 不是疏忽，而是多起点机制的一部分。若要复现具体数值，应改成固定种子并自行重复多次取最优。
>    **例外**：`Code_canonical_model_est_and_fit.m` 中没有 `rng` 调用，在全新 MATLAB 会话中反而是可复现的。
> 4. **bootstrap 非常耗时**（`Readme.txt` 原文：*"the bootstrap takes time"*）。IR 的 bootstrap 还额外降低了迭代与抽样数（`Code_IR_cons_nonparametric_bootstrap.m:6-7`：*"Due to time constraints the number of iterations and draws is lower than in the main estimation"*）——这意味着那几张图的置信带**不能与主估计的精度直接比较**。

**仓库中已发现的不一致（如实记录，未作推测）**：

| 位置 | 问题 |
|---|---|
| `Readme.txt:58` | 提到的 `Code_consumption_est_FE_alternative_normalization.m` **在目录中不存在**；实际文件是 `Code_consumption_est_FE.m`（内容与描述一致） |
| `Code_consumption_parametric_bootstrap.m:39` vs `_fit.m:13-14` | 前者 `Nboot=20`，后者 `Nboot=200; ind=(1:200)`；代码中看不到拼接逻辑，**200 个 replication 的来源不清楚** |
| `Code_consumption_fit_FE.m:6-7` | 代码注释原文：*"In the online version of the supplement, Figure S25b has the wrong scale. The online version of Figure S25c has the correct scale."* |
| `Code_earnings_est.m:396-397` | 归一化只施加于 $\varepsilon$，$\eta$ 的转移系数与初期条件系数**无**归一化行 |
| `Code_earnings_fit.m:13` | `global Vect Vect_dep xx bdw tau` 中的 `xx`/`bdw` 只被未被调用的 `fun_qrlocal.m` 使用，属遗留声明 |
| `floor_par=1` | 代码注释称"=1 in the results in the paper"，但原文 Figure 8 注写的是 $a_{it}\ge 0$。代码的下限加在**对数**资产上（`max(..., floor_par)`），单位与原文表述不一致，复现时需留意 |

**数据可得性**：PSID 样本可由 `data.dta` + `first_stage_regs.do` 完整重建；**挪威数据不可得**（原文附录 C2 描述，2,873 户平衡子样本 2000–2005）。因此 Figure 1b、S16–S19 无法复现。

---

## 7. 关键数值汇总

**样本与网格**

| 量 | 值 |
|---|---|
| 样本 | PSID 工作家庭，1999–2009 双年，平衡面板 |
| $N\times T$ | 792 × 6（= 4,752 观测） |
| 分位点 | $\tau\in\{1/12,\dots,11/12\}$，$L=11$ |
| 年龄范围 | 25–60（均值 44.71，标准差 7.77） |
| 模拟扩样本 | `Nsim=20` → 15,840 个体 |
| 脉冲响应模拟 | 100,000 个体，起始年龄 35，冲击年龄 37 |

**估计参数（`data_hermite2.mat` 实测）**

| 参数 | 值/维数 | 含义 |
|---|---|---|
| `Resqfinal` | 12 × 11 | $\eta_t\mid\eta_{t-1},\text{age}$ 系数 |
| `Resqfinal_e0` | 3 × 11 | $\eta_1\mid\text{age}_1$ 系数 |
| `Resqfinal_eps` | 3 × 11 | $\varepsilon\mid\text{age}$ 系数 |
| `b1, bL` | 4.156, 8.471 | $\eta$ 的 $\lambda^-,\lambda^+$ |
| `b1_e0, bL_e0` | 5.242, 4.047 | $\eta_1$ 的 $\lambda^-,\lambda^+$ |
| `b1_eps, bL_eps` | 2.927, 3.904 | $\varepsilon$ 的 $\lambda^-,\lambda^+$ |

**消费模型结构与自洽性校验**

| 参数 | 值 |
|---|---|
| 消费规则基 | $(a,\eta,\varepsilon,\text{age})$，$(M_1,M_2,M_3,M_4)=(2,2,1,1)$ |
| 资产规则基 | $(a_{t-1},c_{t-1},y_{t-1},\eta_{t-1},\text{age}_t)$，$R=(1,1,1,1,1)$ |
| 期初资产基 | $(\eta_1,\text{age}_1)$，$M_5=M_6=1$ |
| `b1_a, bL_a` | 1.868, 2.044 |

> [!check] 直接读取 `.mat` 确认的维数自洽
> - `Resfinal` = **37×1** = $3\times3\times2\times2$ 个消费规则系数 + 1 个方差 ✔
> - `Resfinal_a2` = **33×1** = $2^5$ 个资产规则系数 + 1 个方差 ✔
> - `Resfinal_a1` = **4×11** = $(M_5{+}1)(M_6{+}1)$ 个系数 × 11 个分位点 ✔
>
> 这从数据侧独立验证了 §5.3 对模型设定的解读，也与论文脚注 35（*"消费规则估计使用次数为 $(2,2,1)$ 的 Hermite 张量积"*）一致。

---

## 8. 学习与迁移建议

1. **最值得迁移的是"用条件分位数刻画 Markov 过程"这一招**。它绕开对转移密度的参数化，允许持续性随状态和冲击大小变化，而 M 步只是分位数回归——计算上极其友好。任何"持续性可能非恒定"的动态面板问题都可考虑。

2. **注意 M 步"不是似然"这一点**。ABB 明确说他们用一串分位数回归替代似然最大化。这带来鲁棒性和计算便利，代价是失去 EM 的单调收敛保证，也因此需要多起点择优（见原文 §5.2 关于"卡在局部状态"的说明）。

3. **消费规则的"对 $\tau$ 可加"设定是被低估的细节**。它使一个分位数模型退化为正态位置—尺度模型，从而 OLS 即为有效估计。好处是简约；代价是**消费分布的形状（偏度、峰度）被假设掉了**——只保留了位置和尺度。若关心消费风险的高阶特征，就需要转向 ABB+Light (2024) 那种估计完整条件分布的做法。这可能是该框架最清晰的边界。

4. **收入侧与消费侧的测量误差处理不对称**：收入侧承认 $\varepsilon$ 与测量误差不可分；消费侧承认 $\nu$ 部分吸收测量误差。两者都是**承认而非解决**——这是使用该框架解读结果时必须记住的。

5. **作为研究起点**：该框架在 CFPS/CHIP 等中国数据上仍有空间——中国的非线性持续性是否存在、房产在资产组合中的角色如何改变保险系数、以及城镇化/流动人口是否产生不同的持续性pattern。

---

## 参考

- Arellano, M., Blundell, R., & Bonhomme, S. (2017). Earnings and Consumption Dynamics: A Nonlinear Panel Data Framework. *Econometrica*, 85(3), 693–734.
- Arellano, M., & Bonhomme, S. (2016). Nonlinear Panel Data Estimation via Quantile Regressions. *Econometrics Journal*, 19(3), C61–C94.
- Blundell, R., Pistaferri, L., & Preston, I. (2008). Consumption Inequality and Partial Insurance. *AER*, 98(5), 1887–1921.
- Celeux, G., & Diebolt, J. (1993). The SEM Algorithm: A Probabilistic Teacher Algorithm Derived from the EM Algorithm for the Mixture Problem.
- Hu, Y., & Schennach, S. (2008). Instrumental Variable Treatment of Nonclassical Measurement Error Models. *Econometrica*, 76(1), 195–216.
- Kim, T.-H., & White, H. (2004). On More Robust Estimation of Skewness and Kurtosis. *Finance Research Letters*, 1(1), 56–73.
- Nielsen, S. F. (2000). The Stochastic EM Algorithm: Estimation and Asymptotic Results. *Bernoulli*, 6(3), 457–489.
- Wei, Y., & Carroll, R. J. (2009). Quantile Regression with Measurement Error. *JASA*, 104(487), 1129–1143.
- 代码仓库内 `Readme.txt` 与 `Simulations/{Calibrated,Estimated}/readme.txt`
