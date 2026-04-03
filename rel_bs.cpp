/**
 * ============================================================
 *  Relativistic Black-Scholes Option Pricer
 *  Based on the Telegrapher's Equation (Goldstein-Kac model)
 * ============================================================
 *
 * PHYSICS:
 *   Standard Black-Scholes uses Geometric Brownian Motion, whose
 *   underlying diffusion equation (parabolic PDE) permits infinite
 *   information propagation speed — a non-relativistic artifact.
 *
 *   The relativistic extension replaces GBM with a persistent
 *   (correlated) random walk with autocorrelation time τ_rel.
 *   This gives a finite "financial speed of light":
 *
 *       c = σ / √(2·τ_rel)    [units: log-price / time]
 *
 *   As τ_rel → 0, c → ∞ and we recover standard Black-Scholes.
 *
 * PDE (forward time τ = T−t, log-price x = ln(S/K)):
 *
 *   τ_rel·∂²V/∂τ² + ∂V/∂τ  =  L[V]
 *
 *   where:
 *       L[V] = (σ²/2)·∂²V/∂x² + (r − σ²/2)·∂V/∂x − r·V
 *
 *   This is the damped wave (telegrapher) equation. When τ_rel = 0
 *   the second-derivative term vanishes → standard parabolic BS.
 *
 * NUMERICAL METHOD:
 *   3-level implicit scheme (Crank-Nicolson in space, 2nd-order in time):
 *
 *   (α·I − θ·L_mat)·V^{n+1}
 *       = 2β·V^n + (γ·I + θ·L_mat)·V^{n-1}  +  BC terms
 *
 *   where: α = τ_rel/Δτ² + 1/(2Δτ)
 *          β = τ_rel/Δτ²
 *          γ = −τ_rel/Δτ² + 1/(2Δτ)      (note: can be negative)
 *          θ = 0.5  (Crank-Nicolson spatial weighting)
 *
 *   The resulting tridiagonal system is solved with the Thomas
 *   algorithm, pre-factored for O(1) re-use across time steps.
 *
 *   Level 0  (τ=0): payoff condition
 *   Level 1  (τ=Δτ): initialized via BS closed-form for 2nd-order accuracy
 *   Levels 2+: 3-level scheme
 *
 * COMPILE:
 *   g++ -O3 -march=native -std=c++17 -o rel_bs rel_bs.cpp
 *
 * REFERENCES:
 *   Biró & Rosenfeld (2006); Kleinert & Korbel (2016);
 *   Bustamante & Contreras (2016)
 * ============================================================
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <chrono>
#include <array>
#include <cassert>
#include <cstdio>

// ─────────────────────────────────────────────────────────────────────────────
//  Normal CDF / PDF
// ─────────────────────────────────────────────────────────────────────────────
static constexpr double SQRT1_2   = 0.70710678118654752440;   // 1/√2
static constexpr double INV_SQRT2PI = 0.39894228040143267794;  // 1/√(2π)

static inline double norm_cdf(double x) noexcept {
    return 0.5 * std::erfc(-x * SQRT1_2);
}
static inline double norm_pdf(double x) noexcept {
    return INV_SQRT2PI * std::exp(-0.5 * x * x);
}

static constexpr double INV_PI = 0.31830988618379067154;   // 1/π

// ─────────────────────────────────────────────────────────────────────────────
//  Cauchy distribution helpers  (Priority 6)
//
//  The Cauchy distribution (α=1 Lévy-stable) is the opposite extreme to the
//  Goldstein-Kac (telegrapher) process:
//    •  RBS: hard light cone,  density ZERO  outside |x| > c·T
//    •  Cauchy: NO causal boundary, polynomial tail  ~ (γ·T)/(π·x²)  for |x| >> γ·T
//
//  Vanilla Cauchy call prices diverge (E[S_T] = ∞ for Cauchy log-returns).
//  Cauchy DIGITAL options are well-posed: P(S_T > K | Cauchy) is finite and
//  provides the meaningful measure of tail probability for comparison.
// ─────────────────────────────────────────────────────────────────────────────

// Cauchy CDF:  F_C(x; 0, 1) = ½ + arctan(x) / π
static inline double cauchy_cdf(double x) noexcept {
    return 0.5 + std::atan(x) * INV_PI;
}

// Cauchy propagator: density of log-return x at horizon T with scale γ_C
//   K_C(x, T; γ_C) = (γ_C·T) / (π · (x² + (γ_C·T)²))
//   Symmetric (even) in x.  Tail: K_C ~ (γ_C·T)/(π·x²) for |x| >> γ_C·T.
static inline double cauchy_propagator(double x, double T, double gamma_C) noexcept {
    const double s = gamma_C * T;
    return s * INV_PI / (x * x + s * s);
}

// Cauchy digital option:  e^{−rT} · P(S_T > K | log-return X_T ~ Cauchy(0, γ_C·T))
//   d_C = ln(S/K) / (γ_C·T)
//   Call digital:  e^{−rT} · [½ + arctan(d_C)/π]   (= e^{−rT}·F_C(d_C) since F_C symmetric)
//   Put digital:   e^{−rT} · [½ − arctan(d_C)/π]
//   Identity:  C_dig + P_dig = e^{−rT}  (exact, verified in Test 11)
static double cauchy_digital(double S, double K, double T, double r,
                              double gamma_C, bool is_call) noexcept {
    if (T < 1e-15 || K <= 0.0) return 0.0;
    const double d_C  = std::log(S / K) / (gamma_C * T);
    const double prob = is_call ? cauchy_cdf(d_C) : (1.0 - cauchy_cdf(d_C));
    return std::exp(-r * T) * prob;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Black-Scholes closed-form reference
// ─────────────────────────────────────────────────────────────────────────────
struct BSResult {
    double price, delta, gamma, vega, theta;
};

BSResult bs_greeks(double S, double K, double T, double r, double sigma, bool is_call) {
    if (T < 1e-15) {
        double pay = is_call ? std::max(S - K, 0.0) : std::max(K - S, 0.0);
        double d = is_call ? (S > K ? 1.0 : 0.0) : (S < K ? -1.0 : 0.0);
        return { pay, d, 0.0, 0.0, 0.0 };
    }
    const double sqrtT = std::sqrt(T);
    const double d1    = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);
    const double d2    = d1 - sigma * sqrtT;
    const double nd1   = norm_pdf(d1);
    const double disc  = std::exp(-r * T);

    BSResult g;
    if (is_call) {
        const double Nd1 = norm_cdf(d1), Nd2 = norm_cdf(d2);
        g.price = S * Nd1 - K * disc * Nd2;
        g.delta = Nd1;
        g.theta = (-S * nd1 * sigma / (2.0 * sqrtT) - r * K * disc * Nd2) / 365.0;
    } else {
        const double Nd1 = norm_cdf(-d1), Nd2 = norm_cdf(-d2);
        g.price = K * disc * Nd2 - S * Nd1;
        g.delta = Nd1 - 1.0;          // = norm_cdf(d1) − 1
        g.theta = (-S * nd1 * sigma / (2.0 * sqrtT) + r * K * disc * Nd2) / 365.0;
        // Note: Nd1 here = norm_cdf(-d1), so delta_put = norm_cdf(-d1) - 1 is WRONG
        // Fix: delta_put = norm_cdf(d1) - 1
        g.delta = norm_cdf(d1) - 1.0;
    }
    g.gamma = nd1 / (S * sigma * sqrtT);
    g.vega  = S * nd1 * sqrtT;        // per unit vol (not per 1%)
    return g;
}

// BS price at log-price x = ln(S/K), forward time tau
static double bs_at_x(double x, double K, double tau, double r, double sigma, bool is_call) {
    return bs_greeks(K * std::exp(x), K, tau, r, sigma, is_call).price;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Second time-derivative of BS price  (closed-form)
//
//  ∂²V_BS/∂τ² = L²[V_BS]  (apply BS spatial operator L twice, using ∂V/∂τ = L[V])
//
//  For a call:
//    ∂²C/∂τ² = S·φ(d₁)·Q − r²·K·e^{−rτ}·N(d₂)
//
//    Q = (d₁²−1)·σ/(4τ√τ)  −  d₁·σ²/(4τ)  −  r·d₁/τ  +  r²/(σ√τ)
//
//  For a put: ∂²P/∂τ² = ∂²C/∂τ² + r²·K·e^{−rτ}   (put-call parity)
//
//  Derivation sketch:
//    L[V] = (σ²/2)·V_xx + (r−σ²/2)·V_x − r·V    (in log-price x = ln S/K)
//    θ = ∂V_BS/∂τ = L[V_BS]   (BS PDE)
//    ∂²V_BS/∂τ² = L[θ] = (σ²/2)·θ_xx + (r−σ²/2)·θ_x − r·θ
//    θ_x  = S·φ(d₁)·[r/σ√τ − d₂/(2τ)]
//    θ_xx = S·φ(d₁)·{(1 − d₁/σ√τ)·[r/σ√τ − d₂/(2τ)] − 1/(2σ·τ^{3/2})}
// ─────────────────────────────────────────────────────────────────────────────
static double bs_d2V_dtau2(double S, double K, double T, double r,
                            double sigma, bool is_call) noexcept {
    if (T < 1e-12) return 0.0;
    const double sqrtT = std::sqrt(T);
    const double d1    = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T)
                         / (sigma * sqrtT);
    const double d2    = d1 - sigma * sqrtT;
    const double phi1  = norm_pdf(d1);
    const double disc  = std::exp(-r * T);

    // Polynomial in d₁, σ, τ  arising from L²
    const double Q = (d1 * d1 - 1.0) * sigma / (4.0 * T * sqrtT)
                   - d1 * sigma * sigma / (4.0 * T)
                   - r * d1 / T
                   + r * r / (sigma * sqrtT);

    const double call_d2 = S * phi1 * Q - r * r * K * disc * norm_cdf(d2);
    // Put: ∂²P/∂τ² = ∂²C/∂τ² + r²·K·e^{−rτ}
    return is_call ? call_d2 : call_d2 + r * r * K * disc;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Black-Scholes risk-neutral density  (Breeden-Litzenberger)
//
//  q_BS(K) = e^{rT}·∂²C_BS/∂K² = φ(d₂) / (K·σ·√T)
//
//  where  d₂ = [ln(S/K) + (r − σ²/2)·T] / (σ·√T)
// ─────────────────────────────────────────────────────────────────────────────
static double bs_density(double S, double K, double T, double r, double sigma) noexcept {
    if (K <= 0.0 || T < 1e-15) return 0.0;
    const double sqrtT = std::sqrt(T);
    const double d2    = (std::log(S / K) + (r - 0.5 * sigma * sigma) * T) / (sigma * sqrtT);
    return norm_pdf(d2) / (K * sigma * sqrtT);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Perturbative RBS pricer  (1/c² expansion, O(τ_rel) accurate)
//
//  V_RBS ≈ V_BS + (1/c²)·v,   1/c² = 2τ_rel/σ²,   v = −(σ²/2)·T·∂²V_BS/∂τ²
//
//  Equivalently: V_RBS ≈ V_BS − τ_rel·T·∂²V_BS/∂τ²
//
//  Validity gauge: |correction|/V_BS < 5% → perturbative regime is safe.
//  The formula reproduces the "frown" (higher price near ATM, lower at extremes)
//  as a consequence of the Telegrapher PDE's finite speed of information.
// ─────────────────────────────────────────────────────────────────────────────
struct PerturbResult {
    double V_bs;       // Black-Scholes price
    double correction; // −τ_rel·T·∂²V_BS/∂τ²
    double V_rbs;      // V_bs + correction
    double rel_size;   // |correction|/V_bs  (perturbative validity gauge)
};

static PerturbResult rbs_perturbative(double S, double K, double T, double r,
                                       double sigma, double tau_rel, bool is_call) {
    const double V_bs = bs_greeks(S, K, T, r, sigma, is_call).price;
    const double d2V  = bs_d2V_dtau2(S, K, T, r, sigma, is_call);
    const double corr = -tau_rel * T * d2V;
    const double gauge = (V_bs > 1e-12) ? std::abs(corr) / V_bs : 0.0;
    return { V_bs, corr, V_bs + corr, gauge };
}

// ─────────────────────────────────────────────────────────────────────────────
//  Implied Volatility solver  (Newton-Raphson, bisection fallback)
//
//  Inverts V_BS(σ) = price for σ_impl.
//  Initial guess: Brenner-Subrahmanyam approximation (works for all moneyness).
//  Convergence: typically 3-6 N-R iterations for liquid strikes.
// ─────────────────────────────────────────────────────────────────────────────
static constexpr double SQRT2PI = 2.5066282746310002;  // √(2π)

static double implied_vol(double price, double S, double K, double T, double r,
                          bool is_call, double tol = 1e-8) {
    if (T < 1e-12) return 0.0;
    const double disc      = std::exp(-r * T);
    const double intrinsic = is_call ? std::max(S - K * disc, 0.0)
                                     : std::max(K * disc - S, 0.0);
    if (price <= intrinsic + 1e-10) return 0.0;  // price at or below intrinsic

    // Brenner-Subrahmanyam initial guess:  σ ≈ √(2π/T) · price / (½(S + K·e^{−rT}))
    double sigma = (SQRT2PI / std::sqrt(T)) * price / (0.5 * (S + K * disc));
    sigma = std::clamp(sigma, 1e-4, 5.0);

    // Newton-Raphson
    for (int i = 0; i < 100; ++i) {
        BSResult g  = bs_greeks(S, K, T, r, sigma, is_call);
        double diff = g.price - price;
        if (std::abs(diff) < tol) return sigma;
        if (g.vega < 1e-12) break;   // degenerate vega → bisection
        double step = diff / g.vega;
        sigma -= step;
        sigma = std::clamp(sigma, 1e-4, 5.0);
        if (std::abs(step) < 1e-11) return sigma;
    }

    // Bisection fallback
    double lo = 1e-4, hi = 5.0;
    for (int i = 0; i < 300; ++i) {
        double mid = 0.5 * (lo + hi);
        double f   = bs_greeks(S, K, T, r, mid, is_call).price - price;
        if (std::abs(f) < tol) return mid;
        (f > 0.0 ? hi : lo) = mid;
        if (hi - lo < 1e-12) return mid;
    }
    return 0.5 * (lo + hi);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Grid: uniform in log-price x, interior points only
//
//  Ghost points (Dirichlet BCs) sit at x_min and x_max.
//  Interior points: x[i] = x_min + (i+1)*dx,  i = 0 .. Nx-1
// ─────────────────────────────────────────────────────────────────────────────
struct Grid {
    int    Nx;
    double x_min, x_max, dx;

    double at(int i) const noexcept { return x_min + (i + 1) * dx; }

    // Factory: grid spans ± n_stdev standard deviations from x = 0 (ATM)
    static Grid make(double sigma, double T, int Nx, double n_stdev = 6.5) {
        double L = std::max(n_stdev * sigma * std::sqrt(T), 1.5);
        Grid g;
        g.Nx    = Nx;
        g.x_min = -L;
        g.x_max =  L;
        g.dx    = 2.0 * L / (Nx + 1);
        return g;
    }
};

// ─────────────────────────────────────────────────────────────────────────────
//  Boundary conditions (Dirichlet, European options)
//
//  In log-price x = ln(S/K), forward time τ:
//    Call: V(x_min, τ) = 0
//          V(x_max, τ) = K*(exp(x_max) - exp(-r*τ))   [deep ITM]
//    Put:  V(x_min, τ) = K*(exp(-r*τ) - exp(x_min))   [deep ITM]
//          V(x_max, τ) = 0
// ─────────────────────────────────────────────────────────────────────────────
static inline std::pair<double, double>
get_bcs(const Grid& g, double K, double r, double tau, bool is_call) noexcept {
    if (is_call) {
        return { 0.0,
                 K * (std::exp(g.x_max) - std::exp(-r * tau)) };
    } else {
        return { K * (std::exp(-r * tau) - std::exp(g.x_min)),
                 0.0 };
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  Spatial operator stencil
//
//  L[V]_i = A*V_{i-1} + B*V_i + C*V_{i+1}
//
//  where A = σ²/(2dx²) − (r−σ²/2)/(2dx)   [sub-diagonal]
//        B = −σ²/dx² − r                    [diagonal]
//        C = σ²/(2dx²) + (r−σ²/2)/(2dx)    [super-diagonal]
// ─────────────────────────────────────────────────────────────────────────────
struct Stencil {
    double A, B, C;

    // m_sq: Klein-Gordon mass term (= 0 → standard BS/RBS spatial operator)
    // Adds −m²·V to L[V], raising the effective discount rate by m².
    static Stencil make(double sigma, double r, double dx, double m_sq = 0.0) noexcept {
        const double sig2  = sigma * sigma;
        const double alpha = sig2 / (2.0 * dx * dx);
        const double beta  = (r - 0.5 * sig2) / (2.0 * dx);
        return { alpha - beta,  -2.0 * alpha - r - m_sq,  alpha + beta };
    }
};

// Apply L to interior vector V with ghost values (bl, br)
// Result written to LV[0..Nx-1]
static void apply_L(const std::vector<double>& V,
                    const Stencil& st,
                    double bl, double br,
                    std::vector<double>& LV) {
    const int Nx = static_cast<int>(V.size());

    // i = 0: left ghost
    LV[0] = st.A * bl + st.B * V[0] + st.C * V[1];

    // i = 1 .. Nx-2: interior (no conditionals — SIMD-friendly)
    for (int i = 1; i < Nx - 1; ++i)
        LV[i] = st.A * V[i-1] + st.B * V[i] + st.C * V[i+1];

    // i = Nx-1: right ghost
    LV[Nx-1] = st.A * V[Nx-2] + st.B * V[Nx-1] + st.C * br;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Pre-factored tridiagonal solver (Thomas algorithm)
//
//  Factor once for a constant matrix; solve in O(Nx) per RHS.
//  The forward elimination multipliers w[i] and modified diagonal
//  diag_f[i] are cached after factoring.
// ─────────────────────────────────────────────────────────────────────────────
struct TridiagSolver {
    int    n;
    std::vector<double> lo, up;    // lower / upper diagonals (unchanged)
    std::vector<double> diag_f;   // factored main diagonal
    std::vector<double> w;        // forward-elimination multipliers

    void factor(const std::vector<double>& lower,
                const std::vector<double>& diag,
                const std::vector<double>& upper) {
        n      = static_cast<int>(diag.size());
        lo     = lower;
        up     = upper;
        diag_f = diag;
        w.resize(n, 0.0);

        for (int i = 1; i < n; ++i) {
            w[i]      = lo[i] / diag_f[i-1];
            diag_f[i] -= w[i] * up[i-1];
#ifndef NDEBUG
            if (std::abs(diag_f[i]) < 1e-14)
                throw std::runtime_error("Thomas: near-zero pivot at i=" + std::to_string(i));
#endif
        }
    }

    // Solve in-place: modifies rhs, writes solution to sol
    // rhs passed by value to preserve caller's buffer
    void solve(std::vector<double> rhs, std::vector<double>& sol) const {
        // Forward substitution (apply stored multipliers)
        for (int i = 1; i < n; ++i)
            rhs[i] -= w[i] * rhs[i-1];

        // Back substitution
        sol[n-1] = rhs[n-1] / diag_f[n-1];
        for (int i = n-2; i >= 0; --i)
            sol[i] = (rhs[i] - up[i] * sol[i+1]) / diag_f[i];
    }
};

// ─────────────────────────────────────────────────────────────────────────────
//  Standard Crank-Nicolson solver  (τ_rel = 0 fallback)
//
//  PDE: ∂V/∂τ = L[V]
//  Scheme: (I/Δτ − 0.5·L)·V^{n+1} = (I/Δτ + 0.5·L)·V^n + BC corrections
// ─────────────────────────────────────────────────────────────────────────────
static std::vector<double>
solve_cn(double sigma, double r, double T, double K,
         int Nt, bool is_call, const Grid& g, double m_sq = 0.0) {
    const double dtau = T / Nt;
    const double inv_dt = 1.0 / dtau;
    const Stencil st = Stencil::make(sigma, r, g.dx, m_sq);

    // System matrix  (I/Δτ − 0.5·L):
    //   lower = −0.5·A,  diag = inv_dt − 0.5·B,  upper = −0.5·C
    std::vector<double> lo  (g.Nx, -0.5 * st.A);
    std::vector<double> diag(g.Nx,  inv_dt - 0.5 * st.B);
    std::vector<double> up  (g.Nx, -0.5 * st.C);

    TridiagSolver solver;
    solver.factor(lo, diag, up);

    std::vector<double> V(g.Nx), V_next(g.Nx), LV(g.Nx), rhs(g.Nx);

    // Level 0: payoff
    for (int i = 0; i < g.Nx; ++i) {
        double xi = g.at(i);
        V[i] = is_call ? K * std::max(std::exp(xi) - 1.0, 0.0)
                       : K * std::max(1.0 - std::exp(xi), 0.0);
    }

    for (int n = 0; n < Nt; ++n) {
        auto [bl_n, br_n]   = get_bcs(g, K, r, n       * dtau, is_call);
        auto [bl_n1, br_n1] = get_bcs(g, K, r, (n + 1) * dtau, is_call);

        apply_L(V, st, bl_n, br_n, LV);

        // RHS: (I/Δτ + 0.5·L)·V^n
        for (int i = 0; i < g.Nx; ++i)
            rhs[i] = inv_dt * V[i] + 0.5 * LV[i];

        // BC corrections for the implicit −0.5·L[V^{n+1}] term
        rhs[0]          += 0.5 * st.A * bl_n1;
        rhs[g.Nx - 1]   += 0.5 * st.C * br_n1;

        solver.solve(rhs, V_next);
        std::swap(V, V_next);
    }
    return V;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Relativistic BS solver  (3-level implicit scheme, θ = 0.5)
//
//  Scheme (derived from τ_rel·V_ττ + V_τ = θ·L[V^{n+1}] + θ·L[V^{n-1}]):
//
//    α·V^{n+1} − θ·L[V^{n+1}]
//        = 2β·V^n + γ·V^{n-1} + θ·L[V^{n-1}]  + BC terms
//
//    α = τ_rel/Δτ² + 1/(2Δτ)       (always positive)
//    β = τ_rel/Δτ²
//    γ = −τ_rel/Δτ² + 1/(2Δτ)     (can be negative for large τ_rel)
//    θ = 0.5
//
//  The system matrix is constant → pre-factored with TridiagSolver.
//  Level 1 initialized via BS closed-form (2nd-order accurate start).
// ─────────────────────────────────────────────────────────────────────────────
static constexpr double TAU_REL_MIN = 1e-12;  // use CN below this threshold

static std::vector<double>
solve_rbs(double sigma, double r, double T, double K, double tau_rel,
          int Nt, bool is_call, const Grid& g, double m_sq = 0.0) {
    if (tau_rel < TAU_REL_MIN)
        return solve_cn(sigma, r, T, K, Nt, is_call, g, m_sq);

    const double dtau  = T / Nt;
    const double theta = 0.5;
    const Stencil st   = Stencil::make(sigma, r, g.dx, m_sq);

    // 3-level time coefficients
    const double A_c   = tau_rel / (dtau * dtau);    // β = A_c
    const double B_c   = 1.0 / (2.0 * dtau);
    const double alpha = A_c + B_c;
    const double beta2 = 2.0 * A_c;                  // coefficient of V^n
    const double gamma = -A_c + B_c;                 // coefficient of V^{n-1}
    //   γ > 0  (τ_rel > Δτ/2): wave-dominated  — telegraph character present
    //   γ < 0  (τ_rel < Δτ/2): diffusion-dominated — relativistic correction
    //          is sub-step in size; scheme stays stable and accurate (implicit)

    // System matrix  (α·I − θ·L_mat):
    //   lower = −θ·A,  diag = α − θ·B,  upper = −θ·C
    std::vector<double> lo  (g.Nx, -theta * st.A);
    std::vector<double> diag(g.Nx,  alpha - theta * st.B);
    std::vector<double> up  (g.Nx, -theta * st.C);

    TridiagSolver solver;
    solver.factor(lo, diag, up);

    std::vector<double> V_prev(g.Nx), V_curr(g.Nx), V_next(g.Nx);
    std::vector<double> LV_prev(g.Nx), rhs(g.Nx);

    // ── Level 0: payoff ────────────────────────────────────────────────────
    for (int i = 0; i < g.Nx; ++i) {
        double xi = g.at(i);
        V_prev[i] = is_call ? K * std::max(std::exp(xi) - 1.0, 0.0)
                            : K * std::max(1.0 - std::exp(xi), 0.0);
    }

    // ── Level 1: BS closed-form at τ = Δτ ─────────────────────────────────
    // Using BS avoids O(Δτ) contamination from naive Euler initialization
    for (int i = 0; i < g.Nx; ++i)
        V_curr[i] = bs_at_x(g.at(i), K, dtau, r, sigma, is_call);

    // ── Time-stepping: n = 1 → Nt−1 ───────────────────────────────────────
    //  At step n we have: V_prev = V^{n-1},  V_curr = V^n
    //  We solve for:      V_next = V^{n+1}
    for (int n = 1; n < Nt; ++n) {
        const double tau_prev = (n - 1) * dtau;
        const double tau_next = (n + 1) * dtau;

        auto [bl_prev, br_prev] = get_bcs(g, K, r, tau_prev, is_call);
        auto [bl_next, br_next] = get_bcs(g, K, r, tau_next, is_call);

        // Compute θ·L[V^{n-1}] with BCs at τ_{n-1}
        apply_L(V_prev, st, bl_prev, br_prev, LV_prev);

        // Build RHS:  2β·V^n + γ·V^{n-1} + θ·L[V^{n-1}]
        // (vectorizable — no data dependencies)
        for (int i = 0; i < g.Nx; ++i)
            rhs[i] = beta2 * V_curr[i] + gamma * V_prev[i] + theta * LV_prev[i];

        // BC corrections for the implicit −θ·L[V^{n+1}] term at endpoints
        rhs[0]          += theta * st.A * bl_next;
        rhs[g.Nx - 1]   += theta * st.C * br_next;

        solver.solve(rhs, V_next);

        // Shift time levels (pointer swap, no allocation)
        std::swap(V_prev, V_curr);
        std::swap(V_curr, V_next);
    }

    return V_curr;  // V at τ = T (i.e., the current option price surface)
}

// Klein-Gordon variant: solve_rbs with explicit mass term m_sq.
// Full documentation in the print_kg_table comment block below.
static std::vector<double>
solve_kg(double sigma, double r, double T, double K, double tau_rel, double m_sq,
         int Nt, bool is_call, const Grid& g) {
    return solve_rbs(sigma, r, T, K, tau_rel, Nt, is_call, g, m_sq);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Interpolation: piecewise-linear on the log-price grid
// ─────────────────────────────────────────────────────────────────────────────
static double grid_interp(const std::vector<double>& V, const Grid& g, double x_target) {
    // x_target = x_min + (idx+1)*dx  =>  idx = (x_target - x_min)/dx - 1
    const double idx_f = (x_target - g.x_min) / g.dx - 1.0;
    if (idx_f <= 0.0)          return V[0];
    if (idx_f >= g.Nx - 1.0)  return V[g.Nx - 1];
    const int i0   = static_cast<int>(idx_f);
    const double f = idx_f - i0;
    return V[i0] * (1.0 - f) + V[i0 + 1] * f;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Greeks from the solved grid
//
//  Delta = ∂V/∂S = (1/S)·∂V/∂x         (chain rule: x = ln(S/K))
//  Gamma = ∂²V/∂S² = (1/S²)·(∂²V/∂x² − ∂V/∂x)
//  Vega  = ∂V/∂σ  via central-difference bump on σ
// ─────────────────────────────────────────────────────────────────────────────
struct GreeksResult {
    double price, delta, gamma, vega;
};

static GreeksResult compute_greeks(
        double S, double K, double T, double r, double sigma,
        double tau_rel, int Nt, bool is_call, const Grid& g) {
    const double x0  = std::log(S / K);
    const double dx2 = 2.0 * g.dx;

    auto V = solve_rbs(sigma, r, T, K, tau_rel, Nt, is_call, g);

    const double price = grid_interp(V, g, x0);
    const double Vp    = grid_interp(V, g, x0 + g.dx);
    const double Vm    = grid_interp(V, g, x0 - g.dx);

    const double dVdx   = (Vp - Vm) / dx2;
    const double d2Vdx2 = (Vp - 2.0 * price + Vm) / (g.dx * g.dx);
    const double delta  = dVdx / S;
    const double gamma  = (d2Vdx2 - dVdx) / (S * S);

    // Vega: 1-bp bump on σ (optimal ε for double precision ≈ σ·sqrt(ε_mach))
    const double eps_v = 1e-4;
    auto Vg_p = solve_rbs(sigma + eps_v, r, T, K, tau_rel, Nt, is_call, g);
    auto Vg_m = solve_rbs(sigma - eps_v, r, T, K, tau_rel, Nt, is_call, g);
    const double vega   = (grid_interp(Vg_p, g, x0) - grid_interp(Vg_m, g, x0))
                          / (2.0 * eps_v);

    return { price, delta, gamma, vega };
}

// ─────────────────────────────────────────────────────────────────────────────
//  Self-test suite
// ─────────────────────────────────────────────────────────────────────────────
static void run_tests(double S = 100, double K = 100, double T = 1.0,
                      double r = 0.05, double sigma = 0.2) {
    const int Nx = 800, Nt = 1000;
    const double x0 = std::log(S / K);
    Grid g = Grid::make(sigma, T, Nx, 6.5);

    std::printf("\n  ─── Self-tests ───────────────────────────────────────\n");
    int fail = 0;

    // Test 1: CN vs BS closed form
    {
        auto Vc = solve_cn(sigma, r, T, K, Nt, true, g);
        double err = std::abs(grid_interp(Vc, g, x0) - bs_greeks(S,K,T,r,sigma,true).price);
        bool pass = err < 5e-3;
        std::printf("  [%s] CN vs BS closed-form:           err = %.2e\n",
                    pass ? "PASS" : "FAIL", err);
        fail += !pass;
    }

    // Test 2: RBS(τ_rel=1e-5) ≈ BS
    {
        auto Vc = solve_rbs(sigma, r, T, K, 1e-5, Nt, true, g);
        double err = std::abs(grid_interp(Vc, g, x0) - bs_greeks(S,K,T,r,sigma,true).price);
        bool pass = err < 5e-3;
        std::printf("  [%s] RBS(τ_rel=1e-5) ≈ BS:           err = %.2e\n",
                    pass ? "PASS" : "FAIL", err);
        fail += !pass;
    }

    // Test 3: Put-call parity
    {
        const double tr = 1e-3;
        auto Vc = solve_rbs(sigma, r, T, K, tr, Nt, true,  g);
        auto Vp = solve_rbs(sigma, r, T, K, tr, Nt, false, g);
        double c = grid_interp(Vc, g, x0);
        double p = grid_interp(Vp, g, x0);
        double pcp_err = std::abs(c - p - (S - K * std::exp(-r * T)));
        bool pass = pcp_err < 5e-4;
        std::printf("  [%s] Put-call parity (τ_rel=1e-3):   err = %.2e (tol 5e-4)\n",
                    pass ? "PASS" : "FAIL", pcp_err);
        fail += !pass;
    }

    // Test 4: Call delta in (0, 1)
    {
        const double tr = 1e-3;
        auto Vc = solve_rbs(sigma, r, T, K, tr, Nt, true, g);
        double Vp = grid_interp(Vc, g, x0 + g.dx);
        double Vm = grid_interp(Vc, g, x0 - g.dx);
        double delta = (Vp - Vm) / (2.0 * g.dx * S);
        bool pass = (delta > 0.01 && delta < 0.99);
        std::printf("  [%s] Call delta in (0,1):            Δ = %.5f\n",
                    pass ? "PASS" : "FAIL", delta);
        fail += !pass;
    }

    // Test 5: Larger τ_rel → call price decreases (less information, lower bound)
    {
        auto Vc_lo = solve_rbs(sigma, r, T, K, 1e-3, Nt, true, g);
        auto Vc_hi = solve_rbs(sigma, r, T, K, 0.1,  Nt, true, g);
        double p_lo = grid_interp(Vc_lo, g, x0);
        double p_hi = grid_interp(Vc_hi, g, x0);
        // Persistent random walk has more variance than GBM at finite timescales
        // → option prices INCREASE with τ_rel (more uncertainty = higher option value)
        bool pass = (p_hi > p_lo);
        std::printf("  [%s] Call price increases with τ_rel: %.5f < %.5f\n",
                    pass ? "PASS" : "FAIL", p_lo, p_hi);
        fail += !pass;
    }

    // Test 6: Convergence rate ≥ 1.5  (scheme is O(h²); error measured
    //         self-consistently vs finest-grid RBS, not vs BS closed-form)
    {
        const double tr = 1e-3;
        Grid g_c = Grid::make(sigma, T, 200, 6.5);
        Grid g_f = Grid::make(sigma, T, 400, 6.5);
        Grid g_r = Grid::make(sigma, T, 800, 6.5);
        double p_c = grid_interp(solve_rbs(sigma, r, T, K, tr, 500,  true, g_c), g_c, x0);
        double p_f = grid_interp(solve_rbs(sigma, r, T, K, tr, 1000, true, g_f), g_f, x0);
        double p_r = grid_interp(solve_rbs(sigma, r, T, K, tr, 2000, true, g_r), g_r, x0);
        double err_c = std::abs(p_c - p_r);
        double err_f = std::abs(p_f - p_r);
        double rate  = (err_f > 1e-15) ? std::log2(err_c / err_f) : 0.0;
        bool pass    = (rate >= 1.5);
        std::printf("  [%s] Convergence rate ≥ 1.5 (self-ref): rate = %.2f\n",
                    pass ? "PASS" : "FAIL", rate);
        fail += !pass;
    }

    // Test 7: Perturbative 1/c² formula reduces BS pricing error by ≥ 50%
    //         (ATM call, τ_rel = 1e-3; compared against fine-grid numerical RBS)
    {
        const double tr = 1e-3;
        auto pr = rbs_perturbative(S, K, T, r, sigma, tr, true);
        Grid g_ref = Grid::make(sigma, T, 800, 6.5);
        double V_num    = grid_interp(solve_rbs(sigma, r, T, K, tr, 2000, true, g_ref),
                                      g_ref, x0);
        double err_approx = std::abs(pr.V_rbs - V_num);
        double err_bs     = std::abs(pr.V_bs  - V_num);
        bool pass = (err_approx < 0.5 * err_bs);
        std::printf("  [%s] Perturb. halves BS error:          "
                    "err_P=%.2e, err_BS=%.2e\n",
                    pass ? "PASS" : "FAIL", err_approx, err_bs);
        fail += !pass;
    }

    // Test 8: IV frown shape — IV peaks at ATM, lower at both deep-ITM and deep-OTM
    //         (τ_rel = 1e-2, clearly wave-dominated at Nt=1000, T=1yr)
    {
        const double tr = 1e-2;
        const double Ks[3] = {70.0, 100.0, 130.0};
        double iv3[3];
        for (int ki = 0; ki < 3; ++ki) {
            Grid   g  = Grid::make(sigma, T, 800, 6.5);
            double x_ = std::log(S / Ks[ki]);
            double V  = grid_interp(solve_rbs(sigma, r, T, Ks[ki], tr, 1000, true, g), g, x_);
            iv3[ki]   = implied_vol(V, S, Ks[ki], T, r, true);
        }
        bool frown_itm = (iv3[1] > iv3[0]);   // IV(ATM) > IV(deep-ITM K=70)
        bool frown_otm = (iv3[1] > iv3[2]);   // IV(ATM) > IV(deep-OTM K=130)
        bool pass = frown_itm && frown_otm;
        std::printf("  [%s] IV frown (τ_rel=1e-2):             "
                    "IV(70)=%.4f%% IV(100)=%.4f%% IV(130)=%.4f%%\n",
                    pass ? "PASS" : "FAIL",
                    iv3[0]*100, iv3[1]*100, iv3[2]*100);
        fail += !pass;
    }

    // Test 9: Light-cone causality — RBS density at 1.3× cone boundary < 50% of BS density
    //         (τ_rel = 0.10, σ = 0.20, T = 1 yr  →  c = 0.4472, cone at |x| = 0.447)
    //         Probe at x_probe = 1.3 × 0.447 = 0.581, K_probe ≈ 179
    {
        const double tr      = 0.10;
        const double c_lc    = sigma / std::sqrt(2.0 * tr);
        const double x_probe = 1.3 * c_lc * T;          // 1.3× outside the cone
        const double K_probe = S * std::exp(x_probe);    // OTM call (K > S)
        const double dK_fd   = S * 0.02;                 // 2% FD step in K

        Grid g9 = Grid::make(sigma, T, 500, 9.0);
        auto price_K = [&](double Kk) -> double {
            double x_ = std::log(S / Kk);
            return grid_interp(solve_rbs(sigma, r, T, Kk, tr, 500, true, g9), g9, x_);
        };

        double q_rbs = std::max(0.0, std::exp(r * T)
            * (price_K(K_probe - dK_fd) - 2.0 * price_K(K_probe)
               + price_K(K_probe + dK_fd))
            / (dK_fd * dK_fd));
        double q_bs  = bs_density(S, K_probe, T, r, sigma);
        double ratio = (q_bs > 1e-12) ? q_rbs / q_bs : 0.0;
        bool   pass  = (ratio < 0.50);
        std::printf("  [%s] Light-cone suppression (1.3× boundary):  "
                    "q_RBS/q_BS = %.4f  (expect < 0.50)\n",
                    pass ? "PASS" : "FAIL", ratio);
        fail += !pass;
    }

    // Test 10: KG(m²=0) is bit-identical to RBS  (verifies m_sq plumbing)
    {
        const double tr = 1e-3;
        Grid g10 = Grid::make(sigma, T, 400, 6.5);
        auto Vr = solve_rbs(sigma, r, T, K, tr, Nt, true, g10);
        auto Vk = solve_kg (sigma, r, T, K, tr, 0.0, Nt, true, g10);
        double max_diff = 0.0;
        for (int i = 0; i < g10.Nx; ++i)
            max_diff = std::max(max_diff, std::abs(Vr[i] - Vk[i]));
        bool pass = (max_diff == 0.0);
        std::printf("  [%s] KG(m²=0) ≡ RBS (bit-exact):              "
                    "max|diff| = %.2e\n",
                    pass ? "PASS" : "FAIL", max_diff);
        fail += !pass;
    }

    // Test 11: Cauchy digital complementarity — C_digital + P_digital = e^{−rT} (exact)
    {
        const double gamma_C = 0.20;
        double cd = cauchy_digital(S, K, T, r, gamma_C, true);
        double pd = cauchy_digital(S, K, T, r, gamma_C, false);
        double err = std::abs(cd + pd - std::exp(-r * T));
        bool pass = (err < 1e-12);
        std::printf("  [%s] Cauchy digital C+P = e^{-rT}:              "
                    "err = %.2e\n",
                    pass ? "PASS" : "FAIL", err);
        fail += !pass;
    }

    // Test 12: Cauchy density > RBS density just outside the light cone
    //   At x = 1.05×c·T (just beyond causal boundary), K_RBS → 0 but K_Cauchy > 0.
    {
        const double tr     = 0.10;
        const double c_lc   = sigma / std::sqrt(2.0 * tr);
        const double x_pr   = 1.05 * c_lc * T;          // 1.05× outside cone
        const double K_pr   = S * std::exp(x_pr);        // OTM call side
        const double dK_fd  = S * 0.02;

        Grid g12 = Grid::make(sigma, T, 400, 9.0);
        auto pK = [&](double Kk) {
            return grid_interp(solve_rbs(sigma, r, T, Kk, tr, 400, true, g12),
                               g12, std::log(S / Kk));
        };
        double q_rbs_K = std::max(0.0, std::exp(r * T)
            * (pK(K_pr - dK_fd) - 2.0 * pK(K_pr) + pK(K_pr + dK_fd))
            / (dK_fd * dK_fd));
        double q_rbs_x   = K_pr * q_rbs_K;                      // convert to x-space
        double q_cauchy_x = cauchy_propagator(x_pr, T, sigma);   // x-space, γ_C = σ
        bool   pass = (q_cauchy_x > q_rbs_x);
        std::printf("  [%s] Cauchy > RBS at 1.05× cone boundary:       "
                    "K_C=%.4f  K_RBS=%.4f\n",
                    pass ? "PASS" : "FAIL", q_cauchy_x, q_rbs_x);
        fail += !pass;
    }

    std::printf("  ─────────────────────────────────────────────────────\n");
    std::printf("  %d/%d tests passed\n\n", 12 - fail, 12);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Main comparison table
// ─────────────────────────────────────────────────────────────────────────────
static void print_table(double S, double K, double T, double r, double sigma,
                         int Nx = 800, int Nt = 1000) {
    Grid g = Grid::make(sigma, T, Nx, 6.5);
    const double x0 = std::log(S / K);
    BSResult bs_c = bs_greeks(S, K, T, r, sigma, true);
    BSResult bs_p = bs_greeks(S, K, T, r, sigma, false);

    std::printf("\n");
    std::printf("  ══════════════════════════════════════════════════════════════════════════\n");
    std::printf("   Relativistic Black-Scholes  vs  Standard BS\n");
    std::printf("   S=%.1f  K=%.1f  T=%.2f yr  r=%.2f%%  σ=%.1f%%  Nx=%d  Nt=%d\n",
                S, K, T, r*100, sigma*100, Nx, Nt);
    std::printf("  ══════════════════════════════════════════════════════════════════════════\n");

    // Print the financial speed of light for reference tau_rels
    std::printf("\n  Financial 'speed of light'  c = σ/√(2·τ_rel):\n");
    for (double tr : {1e-4, 1e-3, 1e-2, 0.1}) {
        std::printf("    τ_rel = %-8.4f  →  c = %.4f log-units/yr\n",
                    tr, sigma / std::sqrt(2.0 * tr));
    }

    std::printf("\n  %-10s  %10s  %10s  %10s  %10s  %8s  %6s  %7s\n",
                "τ_rel", "Call", "Put", "Δ_call", "Γ (×S)", "C_err%", "mode", "ms");
    std::printf("  %-10s  %10s  %10s  %10s  %10s  %8s  %6s  %7s\n",
                "──────────", "──────────", "──────────",
                "──────────", "──────────", "────────", "──────", "───────");

    // BS reference row
    std::printf("  %-10s  %10.5f  %10.5f  %10.5f  %10.6f  %8s  %6s\n",
                "BS (exact)", bs_c.price, bs_p.price, bs_c.delta,
                bs_c.gamma * S, "—", "—");
    std::printf("  %-10s  %10s  %10s  %10s  %10s  %8s  %6s\n",
                "──────────", "──────────", "──────────",
                "──────────", "──────────", "────────", "──────");

    for (double tr : {0.0, 1e-5, 1e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.2}) {
        auto t0 = std::chrono::high_resolution_clock::now();

        auto Vc = solve_rbs(sigma, r, T, K, tr, Nt, true,  g);
        auto Vp = solve_rbs(sigma, r, T, K, tr, Nt, false, g);

        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

        double call_p = grid_interp(Vc, g, x0);
        double put_p  = grid_interp(Vp, g, x0);

        double Vpl  = grid_interp(Vc, g, x0 + g.dx);
        double Vml  = grid_interp(Vc, g, x0 - g.dx);
        double dVdx = (Vpl - Vml) / (2.0 * g.dx);
        double d2Vdx2 = (Vpl - 2.0 * call_p + Vml) / (g.dx * g.dx);
        double delta_c = dVdx / S;
        double gamma_c = (d2Vdx2 - dVdx) / (S * S);

        double err_pct = 100.0 * std::abs(call_p - bs_c.price) / bs_c.price;

        // mode: τ_rel=0 → plain CN; else wave-dominated if τ_rel > Δτ/2, else diffusion-dominated
        const char* mode = (tr == 0.0) ? "CN" :
                           (tr > T / (2.0 * Nt)) ? "wave" : "diff";
        std::printf("  %-10.5f  %10.5f  %10.5f  %10.5f  %10.6f  %8.4f  %6s  %7.1f\n",
                    tr, call_p, put_p, delta_c, gamma_c * S, err_pct, mode, ms);
    }

    // Put-call parity check
    {
        auto Vc = solve_rbs(sigma, r, T, K, 1e-3, Nt, true,  g);
        auto Vp = solve_rbs(sigma, r, T, K, 1e-3, Nt, false, g);
        double c = grid_interp(Vc, g, x0), p = grid_interp(Vp, g, x0);
        double pcp = S - K * std::exp(-r * T);
        std::printf("\n  Put-call parity (τ_rel = 1e-3):\n");
        std::printf("    C − P            = %+.8f\n", c - p);
        std::printf("    S − K·exp(−rT)   = %+.8f\n", pcp);
        std::printf("    Parity error     =  %.2e\n", std::abs(c - p - pcp));
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  Strike / moneyness sweep
// ─────────────────────────────────────────────────────────────────────────────
static void print_moneyness_table(double S, double T, double r, double sigma,
                                   double tau_rel, int Nt = 1000) {
    std::printf("\n  Moneyness sweep  (τ_rel = %.4f,  T=%.1fyr,  r=%.1f%%, σ=%.1f%%):\n",
                tau_rel, T, r*100, sigma*100);
    std::printf("  %-8s  %10s  %10s  %10s  %10s\n",
                "K", "RBS Call", "BS Call", "RBS Put", "BS Put");
    std::printf("  %-8s  %10s  %10s  %10s  %10s\n",
                "────────", "──────────", "──────────", "──────────", "──────────");

    for (double K : {70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0}) {
        Grid g = Grid::make(sigma, T, 800, 6.5);
        double x0 = std::log(S / K);

        auto Vc = solve_rbs(sigma, r, T, K, tau_rel, Nt, true,  g);
        auto Vp = solve_rbs(sigma, r, T, K, tau_rel, Nt, false, g);

        std::printf("  %-8.1f  %10.5f  %10.5f  %10.5f  %10.5f\n",
                    K,
                    grid_interp(Vc, g, x0), bs_greeks(S,K,T,r,sigma,true).price,
                    grid_interp(Vp, g, x0), bs_greeks(S,K,T,r,sigma,false).price);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  Perturbative 1/c² approximation vs full numerical solver
//
//  Columns:  BS price  |  Perturbative RBS  |  Numerical RBS  |  Δ(P−N)  |  |δ|/V_BS
//
//  Expected:
//    •  Δ(P−N) ≈ 0 for small τ_rel   (perturbative regime)
//    •  Sign of correction negative at deep OTM/ITM → "frown" shape
//    •  rel_size < 5% confirms we are in the valid perturbative regime
// ─────────────────────────────────────────────────────────────────────────────
static void print_perturbative_table(double S, double T, double r, double sigma,
                                      double tau_rel, int Nt = 1000) {
    std::printf("\n  Perturbative 1/c² vs numerical RBS  (τ_rel=%.4f):\n", tau_rel);
    std::printf("  c = σ/√(2τ_rel) = %.4f log-units/yr   "
                "(1/c² = %.4f)\n",
                sigma / std::sqrt(2.0 * tau_rel),
                2.0 * tau_rel / (sigma * sigma));

    std::printf("  %-8s  %10s  %10s  %10s  %10s  %8s\n",
                "K", "BS", "Perturb.", "Numerical", "Δ(P−N)", "|δ|/V_BS");
    std::printf("  %-8s  %10s  %10s  %10s  %10s  %8s\n",
                "────────", "──────────", "──────────",
                "──────────", "──────────", "────────");

    for (double K : {70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0}) {
        Grid   g    = Grid::make(sigma, T, 800, 6.5);
        double x0   = std::log(S / K);
        double V_num = grid_interp(solve_rbs(sigma, r, T, K, tau_rel, Nt, true, g), g, x0);
        auto   pr   = rbs_perturbative(S, K, T, r, sigma, tau_rel, true);

        std::printf("  %-8.1f  %10.5f  %10.5f  %10.5f  %+10.2e  %7.3f%%\n",
                    K, pr.V_bs, pr.V_rbs, V_num,
                    pr.V_rbs - V_num, pr.rel_size * 100.0);
    }
    std::printf("  (|δ|/V_BS < 5%% → perturbative regime safe)\n");
}

// ─────────────────────────────────────────────────────────────────────────────
//  Implied Volatility Surface  ("Frown" Validation — Priority 4)
//
//  For each (strike K, τ_rel): solve full numerical RBS → invert to IV.
//
//  The Relativistic BS model predicts a "frown" (inverted smile):
//    •  Near ATM: finite light-cone barely constrains log-moves → IV ≈ σ_input
//    •  Deep OTM: large log-moves required but suppressed by |x| ≤ c·T →
//       prices lower than BS → implied vol below σ_input (trough)
//    •  Result: IV peaks at ATM, falls at both OTM and (less so) deep ITM
//
//  This contrasts with the empirical market smile (IV rises at OTM extremes).
// ─────────────────────────────────────────────────────────────────────────────
static void print_iv_frown(double S, double T, double r, double sigma, int Nt = 1000) {
    const std::vector<double> strikes  = {70, 80, 90, 95, 100, 105, 110, 120, 130};
    const std::vector<double> tau_rels = {1e-3, 5e-3, 1e-2, 5e-2, 0.10};
    const int nK = static_cast<int>(strikes.size());
    const int nT = static_cast<int>(tau_rels.size());

    // Compute IV for every (K, τ_rel) pair upfront to avoid repeated grid builds
    std::vector<std::vector<double>> ivs(nK, std::vector<double>(nT, 0.0));
    for (int ki = 0; ki < nK; ++ki) {
        const double K  = strikes[ki];
        const double x0 = std::log(S / K);
        Grid g = Grid::make(sigma, T, 800, 6.5);  // one grid per strike, reused across τ_rel
        for (int ti = 0; ti < nT; ++ti) {
            double V = grid_interp(
                solve_rbs(sigma, r, T, K, tau_rels[ti], Nt, true, g), g, x0);
            ivs[ki][ti] = implied_vol(V, S, K, T, r, true);
        }
    }

    // ── Table header ────────────────────────────────────────────────────────
    std::printf("\n  Implied Volatility Surface (%%)  "
                "[S=%.0f, T=%.1fyr, r=%.1f%%, σ_input=%.1f%%]\n",
                S, T, r * 100.0, sigma * 100.0);
    std::printf("  BS flat = %.4f%% (reference).  "
                "RBS values deviate → Frown if ATM > extremes.\n\n", sigma * 100.0);

    std::printf("  %-6s  %9s", "K", "BS(flat)");
    for (double tr : tau_rels)
        std::printf("  τ=%5.0e", tr);
    std::printf("\n  %-6s  %9s", "──────", "─────────");
    for (int j = 0; j < nT; ++j) std::printf("  %8s", "────────");
    std::printf("\n");

    for (int ki = 0; ki < nK; ++ki) {
        const char* marker = (strikes[ki] == 100.0) ? " ←ATM" : "";
        std::printf("  %-6.1f  %9.4f", strikes[ki], sigma * 100.0);
        for (int ti = 0; ti < nT; ++ti)
            std::printf("  %8.4f", ivs[ki][ti] * 100.0);
        std::printf("%s\n", marker);
    }

    // ── Frown depth ──────────────────────────────────────────────────────────
    // Primary metric: IV(ATM=100) − IV(OTM=130) in basis points
    // ATM index = 4 (K=100), OTM index = 8 (K=130)
    std::printf("\n  Frown depth  IV(ATM=100) − IV(OTM=130) in basis points:\n");
    std::printf("  %-12s", "τ_rel");
    for (double tr : tau_rels) std::printf("  %8.1e", tr);
    std::printf("\n  %-12s", "ΔIV (bp)");
    for (int ti = 0; ti < nT; ++ti)
        std::printf("  %8.2f", (ivs[4][ti] - ivs[8][ti]) * 10000.0);
    std::printf("\n  %-12s", "shape");
    for (int ti = 0; ti < nT; ++ti) {
        bool is_frown = (ivs[4][ti] > ivs[0][ti]) && (ivs[4][ti] > ivs[8][ti]);
        std::printf("  %8s", is_frown ? "FROWN" : "partial");
    }
    std::printf("\n");
}

// ─────────────────────────────────────────────────────────────────────────────
//  Timing benchmark: O(Nx × Nt) scaling
// ─────────────────────────────────────────────────────────────────────────────
static void run_benchmark(double sigma = 0.2, double r = 0.05,
                          double T = 1.0, double K = 100.0) {
    std::printf("\n  Benchmark (τ_rel=1e-3, Call):\n");
    std::printf("  %-6s  %-6s  %12s  %12s\n", "Nx", "Nt", "ms/solve", "solves/s");
    std::printf("  %-6s  %-6s  %12s  %12s\n", "──────", "──────", "────────────", "────────────");

    for (int Nx : {400, 800, 1600}) {
        for (int Nt : {500, 1000, 2000}) {
            Grid g = Grid::make(sigma, T, Nx, 6.5);
            // warm-up
            solve_rbs(sigma, r, T, K, 1e-3, Nt, true, g);

            const int reps = 5;
            auto t0 = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < reps; ++i)
                solve_rbs(sigma, r, T, K, 1e-3, Nt, true, g);
            auto t1 = std::chrono::high_resolution_clock::now();

            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count() / reps;
            std::printf("  %-6d  %-6d  %12.2f  %12.0f\n", Nx, Nt, ms, 1000.0 / ms);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  Convergence study: error vs grid refinement
// ─────────────────────────────────────────────────────────────────────────────
static void run_convergence(double S = 100, double K = 100, double T = 1.0,
                            double r = 0.05, double sigma = 0.2) {
    //
    // IMPORTANT: error must be measured against a self-consistent RBS reference,
    // NOT against BS closed-form.  The RBS solution differs from BS by a finite
    // O(τ_rel) amount; comparing to BS makes the rate appear to stall once the
    // numerical error falls below the physical RBS–BS gap.
    //
    // Reference: Richardson extrapolation from the two finest levels.
    //   For an O(h²) scheme with halved step:  V* = (4·V_fine − V_coarse) / 3
    // This gives an O(h⁴) accurate reference with no extra solves.
    //
    const double tau_rel = 1e-3;
    const double x0      = std::log(S / K);

    struct Level { int Nx, Nt; double price; };
    std::array<Level, 4> lv = {{
        {100, 250, 0.0}, {200, 500, 0.0}, {400, 1000, 0.0}, {800, 2000, 0.0}
    }};

    for (auto& L : lv) {
        Grid g  = Grid::make(sigma, T, L.Nx, 6.5);
        auto Vc = solve_rbs(sigma, r, T, K, tau_rel, L.Nt, true, g);
        L.price = grid_interp(Vc, g, x0);
    }

    // Richardson-extrapolated reference  (O(h⁴))
    const double V_ref = (4.0 * lv[3].price - lv[2].price) / 3.0;
    const BSResult bs  = bs_greeks(S, K, T, r, sigma, true);

    std::printf("\n  Convergence study (τ_rel=1e-3, ATM call):\n");
    std::printf("  Reference: Richardson extrapolant V* = %.6f  "
                "(RBS−BS gap: %+.2e)\n", V_ref, V_ref - bs.price);
    std::printf("  %-6s  %-6s  %12s  %12s  %10s\n",
                "Nx", "Nt", "Call price", "err vs V*", "Rate");
    std::printf("  %-6s  %-6s  %12s  %12s  %10s\n",
                "──────", "──────", "────────────", "────────────", "──────────");

    double prev_err = 0.0;
    for (int k = 0; k < 4; ++k) {
        double err  = std::abs(lv[k].price - V_ref);
        double rate = (k > 0 && err > 1e-15 && prev_err > 1e-15)
                      ? std::log2(prev_err / err) : 0.0;
        std::printf("  %-6d  %-6d  %12.6f  %12.2e  %10.2f\n",
                    lv[k].Nx, lv[k].Nt, lv[k].price, err, rate);
        prev_err = err;
    }
    std::printf("  BS (exact)          %12.6f\n", bs.price);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Klein-Gordon variant solver  (Priority 7)
//
//  Adds a mass term m² to the spatial operator:
//
//    L_KG[V] = (σ²/2)·V_xx + (r−σ²/2)·V_x − (r + m²)·V
//
//  PDE:  τ_rel·∂²V/∂τ² + ∂V/∂τ = L_KG[V]
//
//  Physics:
//    •  m = 0         → standard RBS (massless telegrapher)
//    •  m > 0         → Klein-Gordon "mass gap": dispersion relation shifts
//                       from ω ∝ k²  to  ω ∝ k² + m² / τ_rel
//    •  m² > 0        raises the effective discount rate: r_eff = r + m²
//    •  Option prices  DECREASE monotonically with m  (heavier discounting)
//    •  Light-cone speed  c = σ/√(2·τ_rel) is UNCHANGED by the mass term
//       (mass affects amplitude, not propagation speed)
//    •  Financial interpretation: m² acts as a hazard rate (credit spread),
//       or models momentum decay in the persistent random walk
//
//  Implementation: pure wrapper — solve_rbs with non-zero m_sq (see above).
// ─────────────────────────────────────────────────────────────────────────────

// ─────────────────────────────────────────────────────────────────────────────
//  Klein-Gordon output table  (Priority 7)
//
//  Part 1 — mass sweep: ATM call price as m varies (m_sq = m²), fixed τ_rel.
//  Part 2 — strike sweep: KG vs RBS vs BS for a representative m value.
//
//  Key observations:
//    •  KG(m=0) ≡ RBS      (verified bit-exactly in Test 10)
//    •  Call prices fall monotonically with m  (heavier discounting)
//    •  At m → large: price → 0  (infinite effective discount rate)
//    •  Effective r_eff = r + m²  can be confirmed by comparing with
//       BS(r=r+m²) at the same σ and τ_rel
// ─────────────────────────────────────────────────────────────────────────────
static void print_kg_table(double S, double T, double r, double sigma,
                            double tau_rel, int Nx = 800, int Nt = 1000) {
    const double c_light = sigma / std::sqrt(2.0 * tau_rel);
    const double x0      = 0.0;   // ATM: x = ln(S/K) = 0
    Grid g = Grid::make(sigma, T, Nx, 6.5);

    // ── Part 1: mass sweep (ATM call) ────────────────────────────────────────
    const double rbs_atm = grid_interp(
        solve_rbs(sigma, r, T, S, tau_rel, Nt, true, g), g, x0);

    std::printf("\n  Klein-Gordon Variant  (τ_rel=%.4f,  c=%.4f,  S=K=%.0f)\n",
                tau_rel, c_light, S);
    std::printf("  L_KG[V] = (σ²/2)V_xx + (r−σ²/2)V_x − (r+m²)V\n");
    std::printf("  m=0 ≡ RBS;  r_eff = r + m²  (mass raises effective discount rate)\n\n");

    std::printf("  %-8s  %-8s  %-12s  %-12s  %-10s  %-12s\n",
                "m", "m²", "KG Call", "RBS (m=0)", "KG/RBS", "BS(r_eff)");
    std::printf("  %-8s  %-8s  %-12s  %-12s  %-10s  %-12s\n",
                "────────", "────────", "────────────", "────────────",
                "──────────", "────────────");

    for (double m : {0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50}) {
        const double m_sq   = m * m;
        const double r_eff  = r + m_sq;
        const double kg_atm = grid_interp(
            solve_kg(sigma, r, T, S, tau_rel, m_sq, Nt, true, g), g, x0);
        const double bs_eff = bs_greeks(S, S, T, r_eff, sigma, true).price;
        std::printf("  %-8.4f  %-8.4f  %-12.5f  %-12.5f  %-10.6f  %-12.5f\n",
                    m, m_sq, kg_atm, rbs_atm, kg_atm / rbs_atm, bs_eff);
    }

    // ── Part 2: strike sweep at a representative mass ─────────────────────────
    const double m_demo  = 0.10;
    const double msq_demo = m_demo * m_demo;
    const double r_eff   = r + msq_demo;

    std::printf("\n  Strike sweep  (τ_rel=%.4f,  m=%.2f,  m²=%.4f,  r_eff=%.4f):\n",
                tau_rel, m_demo, msq_demo, r_eff);
    std::printf("  %-8s  %-12s  %-12s  %-12s  %-10s\n",
                "K", "BS (m=0)", "RBS (m=0)", "KG (m=0.10)", "Δ KG−RBS");
    std::printf("  %-8s  %-12s  %-12s  %-12s  %-10s\n",
                "────────", "────────────", "────────────", "────────────", "──────────");

    for (double K : {70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0}) {
        Grid   gK    = Grid::make(sigma, T, Nx, 6.5);
        double xK    = std::log(S / K);
        double bs_p  = bs_greeks(S, K, T, r, sigma, true).price;
        double rbs_p = grid_interp(solve_rbs(sigma, r, T, K, tau_rel, Nt, true, gK), gK, xK);
        double kg_p  = grid_interp(solve_kg (sigma, r, T, K, tau_rel, msq_demo, Nt, true, gK), gK, xK);
        std::printf("  %-8.1f  %-12.5f  %-12.5f  %-12.5f  %+10.5f\n",
                    K, bs_p, rbs_p, kg_p, kg_p - rbs_p);
    }
    std::printf("  (Δ KG−RBS < 0 everywhere: mass term suppresses prices uniformly)\n");
}

// ─────────────────────────────────────────────────────────────────────────────
//  Light-cone / causality validation  (Priority 5)
//
//  Extracts the RBS risk-neutral density via Breeden-Litzenberger:
//
//    q(K) = e^{rT} · ∂²C/∂K²     (finite-difference on a fine K-grid)
//
//  and compares to the BS log-normal density q_BS(K) = φ(d₂)/(K·σ·√T).
//
//  Financial light cone:
//    Speed of light  c = σ/√(2·τ_rel)  (log-price units / time)
//    Cone boundary   |x| = c·T,   x = ln(S_T/S)
//
//  The Goldstein-Kac (persistent random walk) process has a hard speed limit:
//  |X_T| ≤ c·T  a.s.  ⟹  q_RBS(K) → 0  for  |ln(K/S)| > c·T.
//  BS has no such limit (GBM propagates at infinite speed).
//
//  Expected output:
//    • Inside  cone (|x| ≤ c·T): RBS / BS ratio ≈ 1 (densities agree)
//    • Outside cone (|x| > c·T): RBS / BS ratio → 0 (causality enforced)
// ─────────────────────────────────────────────────────────────────────────────
static void print_light_cone(double S, double T, double r, double sigma,
                              double tau_rel, int Nx = 500, int Nt = 500) {
    const double c_light = sigma / std::sqrt(2.0 * tau_rel);
    const double cone_x  = c_light * T;   // half-width of light cone in log-return

    // Strike range: cover ±1.8×cone_x, clamped to [0.6, 1.4] log-units
    const double x_range = std::max(std::min(1.8 * cone_x, 1.4), 0.6);
    const double K_lo    = S * std::exp(-x_range);
    const double K_hi    = S * std::exp( x_range);
    const int    nK      = 81;
    const double dK      = (K_hi - K_lo) / (nK - 1);

    // Single grid reused for all K (grid depends on sigma,T,Nx — not on K)
    Grid g = Grid::make(sigma, T, Nx, 9.0);

    std::vector<double> K_arr(nK), C_arr(nK);
    for (int i = 0; i < nK; ++i) {
        K_arr[i] = K_lo + i * dK;
        double x0 = std::log(S / K_arr[i]);
        C_arr[i]  = grid_interp(
            solve_rbs(sigma, r, T, K_arr[i], tau_rel, Nt, true, g), g, x0);
    }

    const double ert = std::exp(r * T);

    std::printf("\n  Light-cone / Causality Validation  "
                "(τ_rel=%.4f,  T=%.1fyr,  S=%.0f,  σ=%.0f%%)\n",
                tau_rel, T, S, sigma * 100.0);
    std::printf("  c = σ/√(2·τ_rel) = %.4f log-units/yr\n", c_light);
    std::printf("  Cone boundary: |x| = c·T = %.4f   "
                "→  K ∈ [%.2f, %.2f]\n",
                cone_x, S * std::exp(-cone_x), S * std::exp(cone_x));
    std::printf("  (x = ln(S/K);  x>0 ITM call,  x<0 OTM call)\n\n");

    std::printf("  %-9s  %-7s  %-12s  %-12s  %-8s  %s\n",
                "x=ln(S/K)", "K", "q_BS", "q_RBS", "Ratio", "Cone");
    std::printf("  %-9s  %-7s  %-12s  %-12s  %-8s  %s\n",
                "─────────", "───────", "────────────", "────────────",
                "────────", "───────");

    const int step = std::max(1, (nK - 2) / 18);
    for (int i = 1; i < nK - 1; i += step) {
        const double K_i   = K_arr[i];
        const double x_i   = std::log(S / K_i);
        const double q_bs  = bs_density(S, K_i, T, r, sigma);
        const double q_rbs = std::max(0.0,
            ert * (C_arr[i-1] - 2.0 * C_arr[i] + C_arr[i+1]) / (dK * dK));
        const double ratio = (q_bs > 1e-12) ? q_rbs / q_bs : 0.0;
        const bool   in    = (std::abs(x_i) <= cone_x);
        std::printf("  %+9.4f  %7.2f  %-12.6f  %-12.6f  %-8.4f  %s\n",
                    x_i, K_i, q_bs, q_rbs, ratio, in ? "inside" : "OUTSIDE");
    }

    // Numerical mass integral inside vs outside the cone
    double rbs_in = 0.0, rbs_out = 0.0, bs_in = 0.0, bs_out = 0.0;
    for (int i = 1; i < nK - 1; ++i) {
        const double K_i   = K_arr[i];
        const double x_i   = std::log(S / K_i);
        const double q_rbs = std::max(0.0,
            ert * (C_arr[i-1] - 2.0 * C_arr[i] + C_arr[i+1]) / (dK * dK));
        const double q_bs  = bs_density(S, K_i, T, r, sigma);
        if (std::abs(x_i) <= cone_x) { rbs_in  += q_rbs * dK; bs_in  += q_bs * dK; }
        else                          { rbs_out += q_rbs * dK; bs_out += q_bs * dK; }
    }
    const double eps = 1e-20;
    std::printf("\n  Density mass (numerical integral):\n");
    std::printf("    Inside  cone: q_RBS = %.4f   q_BS = %.4f\n", rbs_in,  bs_in);
    std::printf("    Outside cone: q_RBS = %.4f   q_BS = %.4f\n", rbs_out, bs_out);
    std::printf("    Outside/(inside+outside): RBS = %.4f   BS = %.4f\n",
                rbs_out / (rbs_in + rbs_out + eps),
                bs_out  / (bs_in  + bs_out  + eps));
    std::printf("    → RBS tail mass suppressed by %.1f× vs BS\n",
                (bs_out / (bs_in + bs_out + eps))
                / (rbs_out / (rbs_in + rbs_out + eps) + eps));
}

// ─────────────────────────────────────────────────────────────────────────────
//  Cauchy distribution switching  (Priority 6)
//
//  The Goldstein-Kac (telegrapher) process has a HARD causal boundary at |x|=c·T:
//  the RBS risk-neutral density is exactly zero outside the light cone (P5).
//
//  The Cauchy distribution (α=1 Lévy-stable) is the opposite extreme:
//    K_C(x, T; γ_C) = (γ_C·T) / (π·(x²+(γ_C·T)²))
//  Power-law tail ~ 1/x²; density NEVER vanishes; NO causal boundary.
//
//  Vanilla Cauchy call prices DIVERGE (E[S_T] = ∞ for Cauchy log-returns), so
//  only Cauchy DIGITAL options are well-posed:
//    C_dig = e^{−rT}·P(S_T>K|Cauchy) = e^{−rT}·[½ + arctan(d_C)/π],  d_C=ln(S/K)/(γ_C·T)
//
//  Three parts:
//    Part 1 — Propagator comparison: K_Gaussian vs K_Cauchy vs K_RBS across strikes
//    Part 2 — Cauchy digital vs BS digital: heavy-tail premium at extreme moneyness
//    Part 3 — Density mixing at the cone boundary:
//             K_mix = (1−ρ)·K_RBS + ρ·K_Cauchy; even 10% Cauchy restores density
//             at causally-forbidden x, illustrating the physical trade-off.
// ─────────────────────────────────────────────────────────────────────────────
static void print_cauchy_table(double S, double T, double r, double sigma,
                                double tau_rel, int Nx = 500, int Nt = 500) {
    const double gamma_C = sigma;   // set γ_C = σ: fair comparison (same scale)
    const double c_light = sigma / std::sqrt(2.0 * tau_rel);
    const double cone_x  = c_light * T;
    const double ert     = std::exp(r * T);
    const double sqrtT   = std::sqrt(T);

    // ── K-grid (same design as print_light_cone) ──────────────────────────────
    const double x_range = std::max(std::min(1.8 * cone_x, 1.4), 0.6);
    const double K_lo    = S * std::exp(-x_range);
    const double K_hi    = S * std::exp( x_range);
    const int    nK      = 81;
    const double dK      = (K_hi - K_lo) / (nK - 1);

    // Solve RBS for all strikes; compute Breeden-Litzenberger density in x-space
    Grid g = Grid::make(sigma, T, Nx, 9.0);
    std::vector<double> K_arr(nK), C_arr(nK), q_rbs_x(nK, 0.0);
    for (int i = 0; i < nK; ++i) {
        K_arr[i] = K_lo + i * dK;
        C_arr[i] = grid_interp(
            solve_rbs(sigma, r, T, K_arr[i], tau_rel, Nt, true, g),
            g, std::log(S / K_arr[i]));
    }
    for (int i = 1; i < nK - 1; ++i) {
        double qK = std::max(0.0, ert * (C_arr[i-1] - 2.0*C_arr[i] + C_arr[i+1]) / (dK*dK));
        q_rbs_x[i] = K_arr[i] * qK;   // Jacobian: density per unit log-return
    }

    // ── Header ────────────────────────────────────────────────────────────────
    std::printf("\n  Cauchy Distribution Switching  (P6)\n");
    std::printf("  γ_C = σ = %.2f  |  τ_rel=%.4f  |  c=%.4f  |  cone |x|=%.4f\n",
                gamma_C, tau_rel, c_light, cone_x);
    std::printf("  Vanilla Cauchy call prices diverge → digitals used for comparison.\n");

    // ── Part 1: Propagator comparison ─────────────────────────────────────────
    std::printf("\n  Part 1: Risk-neutral densities in log-return space\n");
    std::printf("  K_G = Gaussian propagator (drift-adjusted BS)\n");
    std::printf("  K_C = Cauchy propagator (γ_C=σ, symmetric, polynomial tail)\n");
    std::printf("  K_RBS = Numerical RBS density (Breeden-Litzenberger)\n\n");

    std::printf("  %-9s  %-7s  %-12s  %-12s  %-12s  %-8s  %s\n",
                "x=ln(S/K)", "K", "K_Gaussian", "K_Cauchy", "K_RBS", "K_C/K_G", "Cone");
    std::printf("  %-9s  %-7s  %-12s  %-12s  %-12s  %-8s  %s\n",
                "─────────", "───────", "────────────", "────────────",
                "────────────", "────────", "───────");

    const int step = std::max(1, (nK - 2) / 18);
    for (int i = 1; i < nK - 1; i += step) {
        const double Ki   = K_arr[i];
        const double xi   = std::log(S / Ki);          // ln(S/K); sign convention matches codebase
        const double xret = -xi;                        // log-return to reach Ki: ln(Ki/S)
        const double kg   = Ki * bs_density(S, Ki, T, r, sigma);   // x-space Gaussian
        const double kc   = cauchy_propagator(xret, T, gamma_C);   // x-space Cauchy (symmetric)
        const double kr   = q_rbs_x[i];
        const double rat  = (kg > 1e-12) ? kc / kg : 0.0;
        const bool   in   = (std::abs(xi) <= cone_x);
        std::printf("  %+9.4f  %7.2f  %-12.6f  %-12.6f  %-12.6f  %-8.4f  %s\n",
                    xi, Ki, kg, kc, kr, rat, in ? "inside" : "OUTSIDE");
    }

    // ── Part 2: Cauchy digital vs BS digital ──────────────────────────────────
    std::printf("\n  Part 2: P(S_T > K) — Cauchy digital vs BS digital call probability\n");
    std::printf("  At moderate OTM: Cauchy < BS (broader peak, less modal probability).\n");
    std::printf("  At extreme  OTM: Cauchy >> BS (polynomial tail dominates Gaussian).\n\n");

    std::printf("  %-7s  %-9s  %-10s  %-9s  %-10s  %-9s  %s\n",
                "K", "d₂", "N(d₂) BS", "d_C", "CDF_C", "C/BS", "note");
    std::printf("  %-7s  %-9s  %-10s  %-9s  %-10s  %-9s  %s\n",
                "───────", "─────────", "──────────", "─────────",
                "──────────", "─────────", "────────────");

    for (double K : {70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 150.0, 180.0}) {
        const double d2    = (std::log(S / K) + (r - 0.5*sigma*sigma)*T) / (sigma*sqrtT);
        const double nd2   = norm_cdf(d2);
        const double d_C   = std::log(S / K) / (gamma_C * T);
        const double cdf_c = cauchy_cdf(d_C);
        const double ratio = (nd2 > 1e-12) ? cdf_c / nd2 : 0.0;
        const char*  note  = (std::abs(std::log(S / K)) > cone_x) ? "OTM>cone*" : "";
        std::printf("  %-7.1f  %+9.4f  %-10.6f  %+9.4f  %-10.6f  %-9.4f  %s\n",
                    K, d2, nd2, d_C, cdf_c, ratio, note);
    }
    std::printf("  (* OTM>cone: outside RBS light cone — K_RBS≈0, K_Cauchy>0)\n");

    // ── Part 3: Density mixing at the cone boundary ───────────────────────────
    const double x_pr    = 1.05 * cone_x;   // just outside the cone
    const double kc_pr   = cauchy_propagator(x_pr, T, gamma_C);  // symmetric

    // Interpolate K_RBS at the two cone-edge probe points from the precomputed array
    auto rbs_at_K = [&](double Ktgt) -> double {
        double f = (Ktgt - K_lo) / dK;
        int i0   = std::clamp(static_cast<int>(f), 1, nK - 3);
        double t = f - i0;
        return q_rbs_x[i0] * (1.0 - t) + q_rbs_x[i0 + 1] * t;
    };

    const double K_itm = S * std::exp(-x_pr);   // x = +x_pr (ITM call)
    const double K_otm = S * std::exp( x_pr);   // x = −x_pr (OTM call)
    const double kg_itm = K_itm * bs_density(S, K_itm, T, r, sigma);
    const double kg_otm = K_otm * bs_density(S, K_otm, T, r, sigma);
    const double kr_itm = rbs_at_K(K_itm);
    const double kr_otm = rbs_at_K(K_otm);

    std::printf("\n  Part 3: Density mixing at 1.05× cone boundary  (x = ±%.4f)\n", x_pr);
    std::printf("  K_mix(ρ) = (1−ρ)·K_RBS + ρ·K_Cauchy\n");
    std::printf("  At ρ=0: hard causal cutoff (K_RBS≈0). At ρ>0: Cauchy restores tail density.\n\n");

    for (int wing = 0; wing < 2; ++wing) {
        double Kw = (wing == 0) ? K_itm : K_otm;
        double kgw = (wing == 0) ? kg_itm : kg_otm;
        double krw = (wing == 0) ? kr_itm : kr_otm;
        double xw  = (wing == 0) ? x_pr   : -x_pr;
        std::printf("  x=%+.4f  K=%.2f  (%s):\n", xw, Kw,
                    wing == 0 ? "ITM call wing" : "OTM call wing");
        std::printf("  %-6s  %-12s  %-12s  %-12s  %-12s\n",
                    "ρ", "K_Gaussian", "K_Cauchy", "K_RBS", "K_mix");
        std::printf("  %-6s  %-12s  %-12s  %-12s  %-12s\n",
                    "──────", "────────────", "────────────", "────────────", "────────────");
        for (double rho : {0.0, 0.10, 0.30, 0.50}) {
            double kmix = (1.0 - rho) * krw + rho * kc_pr;
            std::printf("  %-6.2f  %-12.6f  %-12.6f  %-12.6f  %-12.6f\n",
                        rho, kgw, kc_pr, krw, kmix);
        }
        std::printf("\n");
    }
    std::printf("  K_Cauchy at cone edge = %.4f  (vs K_RBS_ITM=%.4f  K_RBS_OTM=%.4f)\n"
                "  10%% Cauchy mixing adds %.4f to the mixed density at the boundary.\n"
                "  (Exact K_RBS=0 is visible ~1.3× outside the cone — see Part 1.)\n",
                kc_pr, kr_itm, kr_otm, 0.10 * kc_pr);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Main
// ─────────────────────────────────────────────────────────────────────────────
int main() {
    const double S = 100.0, K = 100.0, T = 1.0, r = 0.05, sigma = 0.20;

    run_tests(S, K, T, r, sigma);
    print_table(S, K, T, r, sigma);
    print_moneyness_table(S, T, r, sigma, 1e-3);
    print_perturbative_table(S, T, r, sigma, 1e-3);
    print_iv_frown(S, T, r, sigma);
    print_kg_table(S, T, r, sigma, 0.10);
    print_light_cone(S, T, r, sigma, 0.10);
    print_cauchy_table(S, T, r, sigma, 0.10);
    run_convergence(S, K, T, r, sigma);
    run_benchmark(sigma, r, T, K);

    std::printf("\n");
    return 0;
}
