# Validation status — ContinuousTimeDP v4

## What is implemented

- Source files: v2 base restored from `ContinuousTimeDP_v2.zip` + v4 robust policy comparison, obstacle/impulse operator, uniform LCP backend and VI/QVI solvers.
- `hjb_lcp_solve` standardized sign and checked residual for CompEcon (`lcpsolve`) and optional Moll (`LCP`) backends.
- Full stationary HJB/VI/QVI residual at **updated** policies and **updated** intervention values.
- Regression suite `tests/test_HJBkit.m` retained plus `tests/test_HJB_v4.m`.
- Illustrative exact obstacle, nonlocal reset QVI, robust Markov-switching stop models.

## Evidence and limitations

- Source and mathematical form inspected in current build session; independent NumPy/SciPy checks are in `reference_v4_check.py`.
- On 2026-10-10 the repaired source was run in MATLAB with CompEcon CEtools: `test_HJBkit` passed all 13 checks, `test_HJB_v4` passed, and `test_HJB_v4_regressions` passed. No runtime/acceleration claims.
- Success on three stylized examples does not amount to verification of all real Moll D models. Fixed-cost 2-asset intervention requires implementing a globally sound model-specific intervention optimizer and boundary conditions.
- CompEcon `lcpsolve` uses a safeguarded Newton/semismooth transformation; convergence is not guaranteed in all settings. Diagnostics should be inspected.
- Approximation errors due to collocation versus FD can differ greatly for non-smooth free boundaries.

## Required MATLAB tests

```matlab
run_setup;
test_HJBkit;
test_HJB_v4;
test_HJB_v4_regressions;
demo_v4;
```

Optional cross-check: place Yuval Tassa's `LCP.m` on the path and run the `moll` backend comparison in `test_HJB_v4`.

## Moll-specific experimental examples

`demo_moll_D_v4` includes source-parameter-matched adaptations of Moll's `option_simple_LCP.m` and `car.m`. They are **not** asserted as validated reproductions. Moll's option reflecting stencil and buy-side extrapolation in `car.m` are not identical to the generic operator/boundary convention used here.

## Repairs verified on 2026-10-10

- Tensor grids and solver metadata now remain scalar structs.
- Collocation VI/QVI residuals apply the operator to basis coefficients; the discount term uses nodal values.
- Ordinary reflecting collocation uses the same Neumann row replacement for policy evaluation and convergence checks, including corner precedence.
- Human-capital controls optimize each upwind quadrant and its zero-drift boundaries using CompEcon `goldenx`, with endpoint comparisons and state constraints.
- The original 36-by-23 human-capital demo converged in 7 iterations with residual 3.58e-9.
- On the 100-by-40 FD grid, human capital converged in 21 iterations, residual 5.23e-9, normalized maximum value difference from Moll's original numerical solution 5.66e-8.
- On the 50-by-30 spline grid, collocation converged in 18 iterations, nodal residual 1.44e-8. Reoptimization on a 55-by-55 interior grid gave continuous HJB residual 3.14e-6 and normalized FD/collocation value difference 0.1524 percent. This is a screening check, not a proof of accuracy.
- The unrepaired source's experimental option/car FD demos converged; their boundary conventions remain different from Moll's originals. The repairs do not claim a new cross-validation of those models.
- Reflecting collocation with an obstacle remains explicitly unsupported. Skiba collocation and the full Moll model collection have not been revalidated for this repaired version.
