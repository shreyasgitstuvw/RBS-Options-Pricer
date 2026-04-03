# Relativistic Black-Scholes Option Pricer

A single-file C++ implementation of the **Relativistic Black-Scholes (RBS)** model, which replaces geometric Brownian motion with a persistent (correlated) random walk governed by the **Telegrapher's (damped-wave) PDE**. The model enforces a finite *financial speed of light* `c = σ/√(2·τ_rel)`: no price information propagates faster than `c`, creating a hard causal boundary in the risk-neutral density. As `τ_rel → 0` the model recovers standard Black-Scholes exactly.

## Model Physics

The Telegrapher's PDE in forward time τ = T − t, log-price x = ln(S/K):

```
τ_rel · ∂²V/∂τ² + ∂V/∂τ = L[V]
L[V] = (σ²/2)·∂²V/∂x² + (r − σ²/2)·∂V/∂x − r·V
```

- **τ_rel → 0**: recovers standard Black-Scholes
- **Financial speed of light**: `c = σ / √(2·τ_rel)`
- **Light cone**: risk-neutral density = 0 for |x| > c·T (proven numerically and via Breeden-Litzenberger)
- **IV frown**: RBS produces a volatility frown (ATM IV > wings) without skew parameters
- References: Biró & Rosenfeld (2006), Kleinert & Korbel (2016), Bustamante & Contreras (2016)

## Build

```bash
g++ -O3 -march=native -std=c++17 -o rel_bs rel_bs.cpp
```

Requires: C++17, standard library only (no external dependencies).

## Run

```bash
./rel_bs
```

Output is approximately 390 lines covering all analysis sections below.
See `sample_output.txt` for a complete reference run.

## Key Numerical Results (S=K=100, T=1yr, r=5%, σ=20%)

| τ_rel | ATM Call | vs BS |
|-------|----------|-------|
| BS (exact) | 10.45058 | — |
| 1e-3 | 10.45297 | +1.95e-3 |
| 0.05 | 10.55063 | +0.96% |
| 0.10 | 10.65089 | +1.92% |
| 0.20 | 10.85432 | +3.86% |

**Convergence** (τ_rel=1e-3, Richardson self-reference): O(h²), rates 2.08 → 2.11 → 2.00

**Benchmark** (Nx=800, Nt=1000): ~10 ms/solve ≈ 100 solves/sec

## Features

### Solvers
- **`solve_cn`** — standard Crank-Nicolson (τ_rel=0 fallback, pre-factored Thomas algorithm)
- **`solve_rbs`** — 3-level implicit scheme for the Telegrapher's PDE (θ=0.5 in space)
- **`solve_rbs_american`** — American options via projected penalty (clamp V ≥ intrinsic each step)
- **`solve_rbs_barrier`** — Down-and-out barrier options (zero grid nodes below barrier each step)
- **`solve_kg`** — Klein-Gordon variant (adds mass term m² to spatial operator)

### Analytics
- **`rbs_perturbative`** — Fast O(1) correction: V_RBS ≈ V_BS − τ_rel · T · ∂²V_BS/∂τ²
- **`implied_vol`** — Newton-Raphson IV solver (Brenner-Subrahmanyam init, bisection fallback)
- **`compute_greeks`** — Price, Δ, Γ, vega via finite-difference bumps
- **`calibrate_tau_rel`** — Golden-section search to fit τ_rel to an observed IV surface

### Output Tables
| Function | Content |
|---|---|
| `print_table` | τ_rel sweep: call/put/Δ/Γ/mode/ms |
| `print_moneyness_table` | Strike sweep: RBS vs BS |
| `print_perturbative_table` | BS / 1/c² correction / numerical RBS |
| `print_iv_frown` | IV surface across K × τ_rel; frown depth in bp |
| `print_ascii_iv_smile` | ASCII bar-chart IV smile for visual review |
| `print_greeks_surface` | Δ and Γ across K × τ_rel |
| `print_light_cone` | Breeden-Litzenberger density: RBS vs BS, inside/outside cone |
| `print_cauchy_table` | Cauchy propagator vs Gaussian vs RBS; digital options; density mixing |
| `print_kg_table` | KG mass sweep: m² effect on pricing |
| `print_calibration_demo` | Round-trip calibration: synthesize IV surface, recover τ_rel |
| `print_american_table` | American put premium; Merton theorem (call = European) |
| `print_barrier_table` | Down-and-out barrier discount vs barrier level |
| `run_convergence` | Richardson-reference O(h²) convergence study |
| `run_benchmark` | Nx × Nt timing grid |

## Test Suite (16/16 PASS)

| # | Test |
|---|---|
| 1 | CN vs BS closed-form (err < 1e-3) |
| 2 | RBS(τ_rel=1e-5) ≈ BS |
| 3 | Put-call parity |
| 4 | Call delta ∈ (0, 1) |
| 5 | Call price increases with τ_rel |
| 6 | Convergence rate ≥ 1.5 (Richardson self-ref) |
| 7 | Perturbative correction halves BS error |
| 8 | IV frown shape (ATM > wings) |
| 9 | Light-cone suppression (density = 0 at 1.3× boundary) |
| 10 | KG(m²=0) ≡ RBS bit-exact |
| 11 | Cauchy digital complementarity C+P = e^{-rT} (exact) |
| 12 | Cauchy density > RBS density at 1.05× cone boundary |
| 13 | American put ≥ European put |
| 14 | American call = European call (Merton, no dividends) |
| 15 | Barrier call ≤ Vanilla call |
| 16 | Barrier(H→0) = Vanilla (barrier below grid) |
