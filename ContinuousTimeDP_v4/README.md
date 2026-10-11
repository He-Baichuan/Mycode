# ContinuousTimeDP v4 — stationary HJB, robust HJB, VI, QVI

**Purpose:** solve *individual* infinite-horizon continuous-time Bellman equations. **Not** a KFE, stationary distribution, transition dynamics or general-equilibrium solver.

## Numerical algorithms (unchanged)

- **`HJBsolve.m`**: monotone finite-difference CTMC generator (when the stencil satisfies the positive-type condition), **one implicit pseudo-time step** per policy update, NO Howard inner loop. For stopping/impulse problems, the implicit linear equation is replaced with a frozen-obstacle LCP.
- **`scsolve_jump.m`**: CompEcon polynomial/spline/linear **collocation + policy evaluation**. For stopping/impulse problems the collocation operator is transformed to nodal values before an LCP solve. High-order collocation is **not guaranteed monotone**.

Both use the same stage optimization, Markov regime switching, finite-activity Poisson jump process, optional robust min-over-nature, and optional intervention operator.

## Installation

Requires **MATLAB R2023b** or comparably modern MATLAB, plus Miranda–Fackler CompEcon `CEtools` on your MATLAB path, including `fundefn`, `funnode`, `funbas`, `gridmake`, `lcpsolve`, `goldenx`.

```matlab
addpath(genpath('D:\MATLAB\CompEcon'));
cd('D:\MATLAB\ContinuousTimeDP_v4');
run_setup;
which HJBsolve -all
which scsolve_jump -all
which lcpsolve -all

test_HJBkit       % inherited v2 tests
test_HJB_v4       % v4 analytic examples
test_HJB_v4_regressions % tensor grids, Neumann residuals, constrained controls
demo_v4           % plots
```

Do not place competing solver versions (v2/v3/v4) earlier in MATLAB's search path.

## Four equation classes

**Ordinary HJB**: `rho*V = max_u {f^u + A^u V}`.

**Robust HJB**: `rho*V = max_u min_h {f^(u,h) + A^(u,h) V}`.

**Stopping VI**: `max{G-V, max_u min_h [f + A V - rho V]} = 0`.

**Impulse-control QVI**: `max{M V-V, max_u min_h [f + A V-rho V]} = 0`.

The impulse operator is `M V(s,j) = max_{xi in Xi(s,j)} [R(s,j,xi)+V(T(s,j,xi),to(s,j,xi))]`.
The supremum is computed **over supplied candidates**, not magically over an unbounded continuous intervention domain. The caller must supply a global target optimizer or candidate coverage. Model-specific restrictions must rule out costless instantaneous intervention loops.

## Unified model interface

Required:

```matlab
model.nregimes=2;
model.rho=0.05;
model.boundary='state_constraint';  % or 'reflecting'
model.candidates=@(S,D,W,j,params) ...;  % struct array
```

Each controller candidate supplies `.x` (N-by-p), `.f` (N-vector), `.g` (N-by-d), optional `.sigma` (N-by-d-by-q), `.feasible` (boolean N-vector), `.q` (switching rate row), and `.jump` (finite-activity jump channels). `D` provides `D.forward`, `D.backward`, `D.central`, `D.grad`, `D.value`. To solve difficult nonconvex/branching problems, supply **all** relevant local and global Hamiltonian candidates: the solvers choose the maximum based on the *discretized* operator. Not all continuous controls can be optimized by an automatic generic FOC.

Optional:

```matlab
model.Q=[-0.02 0.02; 0.03 -0.03];   % fixed Markov generator
model.jump=@jump_callback;          % finite-intensity Poisson jump
model.adversary=@nature_callback;   % max_u min_h, rowwise
model.stop=@stop_reward;            % exogenous stopping obstacle
model.intervention=@impulse;        % endogenous impulse/QVI obstacle
```

`model.stop(S,j,params{:})` returns a scalar or N-vector; alternatively specify a scalar, R-vector, or N-by-R numeric matrix. Use `-Inf` in a regime or node if stopping is forbidden.

`model.adversary(S,D,W,j,controller_candidate,params{:})` returns struct alternatives. Each nature alternative can **override** `.f`, `.g`, `.sigma`, `.q`, `.jump`, `.feasible`, and optionally provide `.x` for the nature distortion. The controller keeps its own `.x`. The stage uses `max_controller min_nature` of **discrete** Hamiltonians. A finite candidate set guarantees robustness only against supplied nature candidates. For continuously varying h, calculate its optimizer in this callback.

`model.intervention(S,j,W,evalfun,params{:})` returns a struct array of **intervention candidates** with:

- `.target`: N-by-d continuous post-intervention states (can differ for every row)
- `.to`: scalar or N-vector destination regime index (default current regime)
- `.reward`: scalar or N-vector immediate net payoff (negative for a cost)
- `.feasible`: optional logical N-vector of economically feasible interventions
- `.action`: optional N-by-p optimizer's chosen intervention control

The shared `evalfun(T,to)` evaluates the **current** value function at arbitrary targets using grid interpolation (FD) or basis evaluation (collocation). Instead of specifying `.target` and `.reward`, a candidate may supply `.value` if it explicitly performs global target optimization using `evalfun`—in that case report corresponding target and action as well for useful outputs. An empty struct array / `[]` means no possible intervention in that regime.

Intervention candidates can modify *both* continuous and discrete state, unlike exogenous Poisson jumps. The economic agent decides **whether** to intervene through the QVI.

## LCP backends and sign conventions

The shared adapter `hjb_lcp_solve.m` accepts **standard sign**:

```
l <= V <= u;
F = B*V - b;
V=l => F>=0; l<V<u => F=0; V=u => F<=0.
```

Its backends are:

```matlab
opts.lcp_backend = 'compecon'; % default, lcpsolve(-B,+b,l,u,V0)
opts.lcp_backend = 'moll';     % if Yuval Tassa LCP.m is installed, LCP(B,-b,l,u,V0,false)
```

These two libraries use **opposite** sign conventions. Never pass `(B,-b)` directly to CompEcon's `lcpsolve`. Verify `info.lcp.residual` and `info.residual_max` rather than trusting backend iteration termination.

## Algorithms for QVI

**FD** one policy update, one implicit LCP per iteration:

```
stage_n = max_u min_h (f^u,h + A^u,h * Vn)
Gn = max( Gstop, M(Vn) )
Bn = diag(rho+1/Delta) - A_n
bn = fn + Vn/Delta
Vn1 = LCP(Bn,bn,Gn)                % frozen intervention value
```

**Collocation** one exact *fixed-policy* evaluation LCP per iteration:

```
V_nodes = Phi*c
Bc = diag(rho)*Phi - Acoeff
Bvalue = Bc / Phi
V_nodes_new = LCP(Bvalue, f, max(Gstop, M(Vnodes)))
```

At each iteration a fresh **nonlinear VI/QVI** residual is evaluated:

```
R = max( obstacle(V)-V, f_opt + A_opt*V-rho*V )
```

The continuation row selection is recomputed at the new V. Value changes **and** this residual must be small before success. A QVI problem may have nonunique or nonconvergent naive fixed-point iterations; `qvi_damping` / `collocation_damping` may help but do not guarantee convergence.

## Returned diagnostics

Legacy output APIs are maintained:

```matlab
[V,x,resid,info,A] = HJBsolve(model,space,nodes,V0,opts);
[cv,nodes,V,x,resid,info] = scsolve_jump(model,space,nodes,V0,[],opts);
```

With obstacles, `resid` is the **full VI/QVI residual** (not just continuation equation). For ordinary collocation with reflecting boundaries, boundary rows report the imposed Neumann residual and interior rows report the HJB residual. `info` is a scalar struct, including on tensor grids, and includes `stop_mask`, `intervention_mask`, `obstacle`, `target`, `destination`, `intervention_choice`, `adversary`, `lcp`, `converged` and `residual_history`.

## Verification examples

- `example_model_optstop_v4`: exact `V(a)=max(10,20a)`, stop if `a<.5`.
- `example_model_qvi_v4`: intervention resets state `a→1-a`, switches regime 1→2 with immediate cost 3; exact constant values `V1=17, V2=20`.
- `example_model_robust_stop_v4`: nature minimizes over Markov switching rates `.1,.2`, regime-1 obstacle 15; exact `V=(15,50/3)`.
- The entire pre-existing v2 `test_HJBkit` suite is retained for regression.

### Known limits (do not claim universal HJB-QVI support)

1. **Stationary infinite horizon only**. No lifecycle terminal-value HJB or transition PDE.
2. **Finite-activity jumps**; not arbitrary infinite-activity Lévy laws.
3. Correlated multidimensional FD diffusion is accepted only if the *positive-type stencil* condition is satisfied. General reflected mixed diffusion and boundary corners require more work.
4. **Collocation + reflecting Neumann + obstacle** is explicitly rejected (needs augmented boundary/active-set equations). The FD solver supports the benchmark reflecting diffusion generator but not all stochastic reflection geometries.
5. Collocation VI/QVI may be ill-conditioned, nonmonotone, dense, and unreliable around kinks. Use FD as benchmark and refine grids.
6. Nonconvex impulse optimization must be handled by adequate intervention candidates or optimizer callbacks; the code does **not** automatically find globally optimal targets.
7. No automatic proof of admissible no-Zeno cost structure; QVI convergence depends on model.
8. The CompEcon and Moll LCP implementations can fail on difficult/large problems; residual flags allow detection, not magic recovery.
9. This is an **HJB-only library**, not an equilibrium or KFE solver.

## Validation status

See `VALIDATION_STATUS.md`. The repaired source passed `test_HJBkit` (13 checks), `test_HJB_v4`, and `test_HJB_v4_regressions` in MATLAB on 2026-10-10. The human-capital FD and collocation models also converged and were checked against Moll's reference and between collocation nodes. These checks do not establish universal convergence for other models.

### Numerical detail: sparse LCP transformation

`hjb_lcp_solve` temporarily configures CompEcon's `lcpsolve` to use the supplied `hjb_lcp_minmax.m` transformation. This reproduces CompEcon's `minmax.m` sign convention but constructs the active-set Jacobian with `speye` rather than a dense identity. The original global CompEcon `lcpsolve` transform setting is restored through `onCleanup`. It avoids potential `0/0` derivatives at exact obstacle kinks in the original `smooth` transformation. The original CompEcon `lcpsolve` Newton iterations remain unchanged.

## Moll original D examples (experimental adapters)

The package additionally contains **parameter-matched models** derived from the official Moll MATLAB scripts:

- `example_model_moll_option_v4.m` + `demo_moll_D_v4.m`: `[0.1,1]` state, `rho=.05`, drift `-.01`, sigma `.01*x`, reward flow `sqrt(x)`, stop reward `10`, reflecting diffusion.
- `example_model_moll_car_v4.m` + `demo_moll_D_v4.m`: two durable ownership states; CRRA=2, interest `.045`, income `.1`, durable flow `.25`, buying asset cost `.2`, selling gain `.1`, asset interval `[-.02,3]`.

These are **real-model interface adapters**, but have **not** been numerically validated in MATLAB against the original Moll output. In particular the Moll option code and this toolkit use slightly different reflecting-boundary stencils, and the Moll car code's boundary extrapolation at forbidden buys is here replaced by a strict feasible-action mask. Compare stopping thresholds, continuation residuals and grid-convergence, not raw operator equality. `demo_moll_D_v4` is not part of the mandatory passing test suite yet.
