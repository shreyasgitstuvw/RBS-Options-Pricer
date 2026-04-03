# Lovable Build Brief — RBS Options Pricer Dashboard

## What to build

A dark-themed financial analytics dashboard called **"RBS Options Pricer"** with a sidebar navigation and 6 panels. Each panel shows:
1. An interactive Recharts chart
2. A styled "Analysis" card with a title and multi-paragraph inference text

The app fetches all data from a local FastAPI server at `http://localhost:8000`.

---

## API

**Base URL:** `http://localhost:8000`

**Health check:** `GET /` → `{ "status": "ok" }`

**Model parameters:** `GET /api/meta` → `{ "params": { "S": 100, "K": 100, "T": 1.0, "r": 0.05, "sigma": 0.20 }, "modes": [...] }`

**Data endpoint:** `GET /api/csv/{mode}` → returns:
```json
{
  "mode": "iv_surface",
  "rows": 130,
  "columns": ["K", "tau_rel", "iv_pct"],
  "data": [
    { "K": 70.0, "tau_rel": 0.0001, "iv_pct": 20.002674 },
    ...
  ],
  "annotations": {
    "title": "Implied Volatility Smile — The RBS Frown",
    "body": "...",
    "key_insight": "..."
  }
}
```

---

## Layout

- **Left sidebar** (fixed, dark `#111`): logo "RBS Pricer", 6 nav items (one per mode), model parameters shown at bottom (S=100, K=100, T=1yr, r=5%, σ=20%)
- **Main content** (scrollable, dark `#0f0f0f`): selected panel fills the right side
- Each panel:
  - Header: mode title from `annotations.title`
  - Chart area: ~500px tall
  - **Key Insight** badge: yellow/amber pill with `annotations.key_insight` text
  - **Analysis** card: dark `#1a1a1a` card, `annotations.body` rendered with paragraph breaks

Use Tailwind. Font: Inter or system-ui. Accent colour: `#4d96ff` (blue). Loading spinner while fetching.

---

## Panel Specifications

### 1. `iv_surface` — IV Smile
**Columns:** `K` (strike), `tau_rel` (relaxation time), `iv_pct` (implied vol %)

**Chart:** `LineChart` from Recharts
- X-axis: `K` (strikes 70–130)
- Y-axis: `iv_pct` (range ~19–21%)
- One `Line` per unique `tau_rel` value (10 lines total)
- Color each line from a blue→purple gradient
- Reference line at y=20 (dashed grey, label "BS flat")
- Legend shows τ values in scientific notation
- Tooltip shows K, IV%, τ_rel

### 2. `density` — Risk-Neutral Density
**Columns:** `x` (log-return), `K`, `k_rbs`, `k_gauss`, `k_cauchy`, `inside_cone` (0 or 1)

**Chart:** `ComposedChart` from Recharts
- X-axis: `x` (log-return, range approx -0.75 to +0.75)
- Y-axis: density value
- Three `Line` components: `k_rbs` (blue, strokeWidth=2), `k_gauss` (green dashed), `k_cauchy` (red dotted)
- `ReferenceArea` for x values where `inside_cone === 0` (fill red, opacity 0.08, label "Forbidden zone")
- Two `ReferenceLine` at the min and max x where `inside_cone` changes (the cone boundaries) — yellow dashed
- Legend: "RBS (Breeden-Litzenberger)", "Gaussian (BS)", "Cauchy (γ_C=σ)"
- Tooltip shows x, all three densities

### 3. `greeks` — Delta and Gamma
**Columns:** `K` (strike), `tau_rel`, `delta`, `gamma`

**Chart:** Two side-by-side `LineChart` components (use CSS grid, 50/50)
- **Left chart — Delta:**
  - X-axis: `K`
  - Y-axis: `delta` (range 0–1)
  - One line per unique `tau_rel` (7 lines), blue→purple gradient
  - Reference line at delta=0.5 (ATM boundary), dashed grey
- **Right chart — Gamma:**
  - X-axis: `K`
  - Y-axis: `gamma`
  - One line per unique `tau_rel` (7 lines), green→yellow gradient
  - Tooltip shows K, delta or gamma, τ_rel

### 4. `american` — American Put Premium
**Columns:** `K`, `tau_rel`, `european_put`, `american_put`, `premium`, `premium_pct`

**Chart:** Two side-by-side `LineChart` components
- **Left chart — Premium ($):**
  - X-axis: `K`
  - Y-axis: `premium`
  - One line per `tau_rel` (4 lines), red→orange gradient
  - Tooltip shows K, premium $, τ_rel
- **Right chart — Premium (%):**
  - X-axis: `K`
  - Y-axis: `premium_pct`
  - One line per `tau_rel` (4 lines), same gradient
  - Y-axis label: "% of European price"
  - Tooltip shows K, premium %, τ_rel

### 5. `barrier` — Barrier Discount
**Columns:** `H` (barrier level), `tau_rel`, `vanilla`, `barrier_price`, `discount`, `discount_pct`

**Chart:** Two side-by-side `LineChart` components
- **Left chart — Price:**
  - X-axis: `H` (barrier level 50–98)
  - Y-axis: `barrier_price`
  - One line per `tau_rel` (4 lines), purple→blue gradient
  - `ReferenceLine` at y = vanilla price (dashed grey, label "Vanilla")
  - Tooltip shows H, barrier price, τ_rel
- **Right chart — Discount %:**
  - X-axis: `H`
  - Y-axis: `discount_pct`
  - One line per `tau_rel` (4 lines)
  - Tooltip shows H, discount %, τ_rel

### 6. `calibration` — RMSE Landscape
**Columns:** `tau_rel`, `rmse_bp`

**Chart:** Single `LineChart` (full width)
- X-axis: `tau_rel` — use log scale. Since Recharts doesn't support log X natively, pre-process: map tau_rel to `log10(tau_rel)` and label ticks as `10^x` (e.g. -5, -4, -3, -2, -1)
- Y-axis: `rmse_bp` (basis points, range 0–25)
- `Line`: blue, strokeWidth=2, no dots (too many points)
- `ReferenceLine` at `log10(0.05) ≈ -1.301`, yellow dashed, label "τ_true = 0.05"
- `ReferenceArea` for x range where `rmse_bp < 1` — green fill, opacity 0.15, label "< 1 bp"
- Tooltip shows τ_rel (formatted as scientific notation), RMSE in bp

---

## Shared Components

### `<AnnotationCard>` component
```
Props: { title: string, body: string, key_insight: string }

Renders:
  - Panel title (h2, white, text-xl font-semibold)
  - Key insight: amber pill badge (bg-amber-900/40 text-amber-300 border border-amber-700)
    with lightning bolt icon + text
  - Analysis card: dark bg (#1a1a1a), rounded-xl, p-6, mt-4
    - "Analysis" label (text-xs text-gray-500 uppercase tracking-wider mb-2)
    - body text split on \n\n, each paragraph as <p class="text-gray-300 text-sm mb-3">
```

### Loading state
Centered spinner (blue, animate-spin) + "Computing..." text while fetch is in progress.

### Error state
Red alert card: "Failed to load data. Is the server running at localhost:8000?"

---

## Sidebar nav items (in order)
1. IV Surface — `#4d96ff` icon (wave)
2. Density — `#6bcb77` icon (bell curve)
3. Greeks — `#ffd93d` icon (delta symbol)
4. American Options — `#ff6b6b` icon (clock)
5. Barrier Options — `#c77dff` icon (wall)
6. Calibration — `#ff9a3c` icon (target)

Active nav item: left border highlight in the item's colour, slight bg highlight.

---

## Parameter display (sidebar bottom)
Small card showing the model parameters:
```
Model Parameters
S = 100  (spot price)
K = 100  (strike)
T = 1 yr
r = 5%
σ = 20%
```
Fetch from `GET /api/meta` on load.

---

## Tech stack
- Vite + React + TypeScript
- Tailwind CSS
- Recharts (for all charts)
- No external state management — simple `useState` + `useEffect` per panel
- Fetch API (no axios needed)
