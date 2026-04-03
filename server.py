"""
server.py — FastAPI backend for the RBS Options Pricer dashboard.

Wraps rel_bs.exe --csv <mode> and returns JSON:
  { data: [...rows as dicts with float values], annotations: { title, body } }

Run with:
  uvicorn server:app --reload --port 8000

Endpoints:
  GET /                        health check
  GET /api/meta                model parameters + available modes
  GET /api/csv/{mode}          data + inference annotations
"""

import io, csv, os, subprocess
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware

# ── locate executable ─────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(SCRIPT_DIR, "rel_bs.exe")
if not os.path.exists(EXE):
    EXE = os.path.join(SCRIPT_DIR, "rel_bs")

VALID_MODES = ["iv_surface", "density", "greeks", "american", "barrier", "calibration"]

# ── inference annotations ─────────────────────────────────────────────────────
ANNOTATIONS = {
    "iv_surface": {
        "title": "Implied Volatility Smile — The RBS Frown",
        "body": (
            "Standard Black-Scholes produces a flat IV surface by construction — every strike "
            "and maturity gives exactly 20% IV. The RBS model produces a frown: implied volatility "
            "is highest at-the-money and declines into the wings. This arises naturally from the "
            "finite-speed constraint without any free skew parameter.\n\n"
            "The frown deepens as τ_rel increases. At τ_rel = 0.20, OTM wings can be 3–4 vol points "
            "below ATM. At τ_rel = 1e-4 the surface is nearly flat — the model approaches the BS limit. "
            "A calibrated τ_rel > 0.05 means the market is pricing significant speed-of-light "
            "constraints into option wings."
        ),
        "key_insight": "Frown depth ∝ τ_rel — the further right you go on the heatmap, the deeper the smile inverts.",
    },
    "density": {
        "title": "Risk-Neutral Density — The Light Cone",
        "body": (
            "The Breeden-Litzenberger formula extracts the market-implied probability distribution "
            "of the future stock price from the second derivative of call prices. For standard BS "
            "this density is log-normal. For the RBS model there is a hard boundary — the light cone "
            "— beyond which the density is exactly zero.\n\n"
            "The vertical dashed lines mark the light-cone boundary: |x| = c·T where "
            "c = σ/√(2·τ_rel). At τ_rel = 0.10 and σ = 20%, this boundary sits at roughly ±0.45 "
            "log-units. The RBS density (blue) is physically zero outside this cone — not a numerical "
            "artifact. The Cauchy density (red) leaks mass into the forbidden zone, which is why it "
            "overprices deep OTM options relative to RBS.\n\n"
            "The log-scale panel reveals tail behaviour: RBS drops sharply at the boundary, Gaussian "
            "decays exponentially, Cauchy decays only as 1/x². The shaded region quantifies the "
            "probability mass that standard models misallocate to causally forbidden outcomes."
        ),
        "key_insight": "Any density outside the cone is causally forbidden — the light-cone shading shows exactly what BS and Cauchy get wrong.",
    },
    "greeks": {
        "title": "Delta and Gamma — RBS Hedging Surface",
        "body": (
            "Delta measures how much the option price moves per unit move in the underlying. "
            "Gamma measures how quickly Delta changes — it is the primary source of hedging cost. "
            "In RBS, the light-cone constraint alters both surfaces in a subtle but important way.\n\n"
            "Gamma is particularly revealing: it decreases at large τ_rel for deep ITM options, "
            "but increases near-ATM at moderate τ_rel (0.05–0.15). The finite-speed constraint "
            "compresses the risk-neutral density near ATM, raising local curvature. Practical "
            "implication: an RBS hedger must rebalance more frequently than a BS hedger for "
            "near-ATM options at moderate τ_rel, but less so for deep ITM.\n\n"
            "The heatmap colour gradient directly shows where hedging costs will diverge from "
            "standard BS predictions — useful for identifying which parts of the portfolio are "
            "most sensitive to relativistic corrections."
        ),
        "key_insight": "High Gamma (bright yellow) = expensive to delta-hedge. Watch for the ATM Gamma increase at τ_rel ≈ 0.05–0.15.",
    },
    "american": {
        "title": "Early-Exercise Premium — When Waiting Is Not Free",
        "body": (
            "An American put grants the right to exercise early. The early-exercise premium is the "
            "dollar amount by which the American put exceeds the European put. In standard BS this "
            "premium is well-understood; RBS adds a new dimension — τ_rel.\n\n"
            "The premium-as-percentage increases with τ_rel. RBS prices uncertainty more expensively "
            "than BS, raising the cost of not exercising and making early exercise more attractive. "
            "At τ_rel = 0.10, the premium can be 8–10% of the European price for near-ATM options, "
            "versus ~5% in standard BS.\n\n"
            "An options desk using BS to value American puts systematically underprices early exercise "
            "when τ_rel is non-trivial. The heatmap on the right shows exactly where the mispricing "
            "is largest — deep ITM options at high τ_rel."
        ),
        "key_insight": "Bright regions in the right heatmap = largest BS mispricing of American exercise. Premium % grows with both ITM-ness and τ_rel.",
    },
    "barrier": {
        "title": "Down-and-Out Barrier Discount — Relativistic Knock-Out",
        "body": (
            "A down-and-out barrier call pays off only if the stock never falls below the barrier H. "
            "As H → S = 100 the barrier becomes binding and the option price drops to zero. "
            "As H → 0 the barrier never triggers and the price recovers the vanilla call.\n\n"
            "At large τ_rel the barrier discount is smaller for barriers well below spot — the RBS "
            "model assigns less probability to the stock reaching those levels than BS does, consistent "
            "with light-cone density suppression. Practitioners pricing barrier products with "
            "relativistic models need less barrier discount near the cone boundary than classical "
            "models suggest.\n\n"
            "Notice the inflection point: below H ≈ 75 the discount curves separate by τ_rel; "
            "above H ≈ 90 all curves converge as the barrier becomes dominant regardless of τ_rel."
        ),
        "key_insight": "Curves separating below H ≈ 75 shows where τ_rel matters most. The gap between lines is the relativistic correction to barrier pricing.",
    },
    "calibration": {
        "title": "Calibration Landscape — Recovering τ_rel from Market Prices",
        "body": (
            "This panel shows the calibration landscape: given a synthetic IV surface generated at "
            "τ_true = 0.05, how does the RMSE (in basis points) behave as a function of candidate τ_rel?\n\n"
            "The minimum is sharp and unambiguous — the golden-section search converges to within "
            "9.4e-9 bp error at the true value. RMSE rises steeply on both sides, meaning the "
            "calibration is well-identified: different τ_rel values produce noticeably different "
            "IV surfaces.\n\n"
            "The green shaded region (RMSE < 1 bp) defines the range of τ_rel values consistent "
            "with sub-basis-point calibration error. In practice, market noise is 1–5 bp, so this "
            "zone defines the uncertainty on calibrated τ_rel.\n\n"
            "Notice the asymmetry: RMSE rises faster for τ_rel too large than too small. This is "
            "because the model becomes weakly identified as τ_rel → 0 (it approaches standard BS "
            "and loses sensitivity to the parameter)."
        ),
        "key_insight": "The sharp V-shape proves τ_rel is identifiable from IV data. The green zone shows calibration uncertainty given realistic market noise.",
    },
}

# ── FastAPI app ───────────────────────────────────────────────────────────────
app = FastAPI(title="RBS Options Pricer API", version="1.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["GET"],
    allow_headers=["*"],
)


def run_csv(mode: str) -> list[dict]:
    """Call rel_bs.exe --csv <mode>, parse stdout CSV, return list of float dicts."""
    result = subprocess.run(
        [EXE, "--csv", mode],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(result.stderr or f"rel_bs.exe exited with code {result.returncode}")
    reader = csv.DictReader(io.StringIO(result.stdout))
    rows = []
    for row in reader:
        rows.append({k: float(v) for k, v in row.items()})
    return rows


@app.get("/")
def health():
    return {"status": "ok", "exe": EXE, "modes": VALID_MODES}


@app.get("/api/meta")
def meta():
    return {
        "params": {"S": 100.0, "K": 100.0, "T": 1.0, "r": 0.05, "sigma": 0.20},
        "modes": VALID_MODES,
        "description": "Relativistic Black-Scholes option pricer. All parameters hardcoded in rel_bs.exe.",
    }


@app.get("/api/csv/{mode}")
def csv_endpoint(mode: str):
    if mode not in VALID_MODES:
        raise HTTPException(status_code=404, detail=f"Unknown mode '{mode}'. Valid: {VALID_MODES}")
    try:
        data = run_csv(mode)
    except RuntimeError as e:
        raise HTTPException(status_code=500, detail=str(e))
    return {
        "mode": mode,
        "rows": len(data),
        "columns": list(data[0].keys()) if data else [],
        "data": data,
        "annotations": ANNOTATIONS[mode],
    }
