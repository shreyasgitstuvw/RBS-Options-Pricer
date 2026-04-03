# Lovable Build Brief — RBS Options Pricer Dashboard (v2)

## What to build

A dark-themed financial analytics dashboard called **"RBS Options Pricer"** with a sidebar and 7 panels. Each panel shows:
1. An interactive **Plotly.js** chart (3D surface or 2D line — specified per panel)
2. A **Key Insight** badge
3. An **Analysis** card with multi-paragraph inference text

All data comes from a local FastAPI server at `http://localhost:8000`.

**Use Plotly.js for all charts** (not Recharts — Recharts has no 3D support). Import via CDN or npm: `plotly.js-dist`.

---

## API

**Base URL:** `http://localhost:8000`

### Two endpoint families

**Flat rows** (for 2D panels): `GET /api/csv/{mode}`
```json
{
  "mode": "density",
  "rows": 99,
  "columns": ["x", "K", "k_rbs", "k_gauss", "k_cauchy", "inside_cone"],
  "data": [{ "x": -0.68, "K": 50.4, "k_rbs": 0.0, "k_gauss": 0.003, "k_cauchy": 0.125, "inside_cone": 0 }, ...],
  "annotations": { "title": "...", "body": "...", "key_insight": "..." }
}
```

**Plotly-ready grids** (for 3D panels): `GET /api/grid/{mode}`
```json
{
  "mode": "iv_surface",
  "x": [70.0, 75.0, 80.0, ..., 130.0],
  "y": [0.001, 0.005, 0.01, 0.05, 0.10, 0.15, 0.20],
  "x_label": "K",
  "y_label": "tau_rel",
  "surfaces": {
    "iv_pct": [[20.002, 20.003, ...], [19.989, ...], ...]
  },
  "annotations": { "title": "...", "body": "...", "key_insight": "..." }
}
```

The `surfaces` object maps column names to 2D arrays of shape `nY × nX`.
Pass directly to Plotly: `{ type: 'surface', x: data.x, y: data.y, z: data.surfaces.iv_pct }`

**Meta:** `GET /api/meta` → `{ "params": { "S": 100, "K": 100, "T": 1.0, "r": 0.05, "sigma": 0.20 }, "modes": [...], "grid_modes": [...] }`

---

## Layout

- **Left sidebar** (fixed 220px, bg `#111827`): logo "RBS Pricer" (white, bold), 7 nav items, model parameters card at bottom
- **Main content** (bg `#030712`, scrollable): selected panel
- Active nav item: left border `3px solid <accent>`, bg `#1f2937`
- Panel structure:
  ```
  [Panel Title h2]
  [Key Insight pill — amber]
  [Chart — Plotly div, height 520px]
  [Analysis card — dark bg, body text]
  ```

Use **Tailwind CSS**. Font: Inter. Show loading spinner while fetching.

---

## Panel Specifications (7 panels)

---

### Panel 1: `iv_surface` — IV Smile (3D Surface)
**Endpoint:** `GET /api/grid/iv_surface`
**Nav accent:** `#4d96ff`

**Chart:** Single Plotly `surface` trace
```javascript
{
  type: 'surface',
  x: data.x,           // K values (13 strikes)
  y: data.y,           // tau_rel values (10 values)
  z: data.surfaces.iv_pct,
  colorscale: 'Plasma',
  colorbar: { title: 'IV (%)' },
  hovertemplate: 'K: %{x}<br>τ_rel: %{y}<br>IV: %{z:.3f}%<extra></extra>'
}
```
Layout:
```javascript
{
  scene: {
    xaxis: { title: 'Strike K' },
    yaxis: { title: 'τ_rel', type: 'log' },
    zaxis: { title: 'IV (%)' },
    camera: { eye: { x: 1.5, y: -1.5, z: 0.8 } }
  },
  paper_bgcolor: '#030712', plot_bgcolor: '#030712',
  font: { color: '#d1d5db' }, margin: { t: 20 }
}
```

---

### Panel 2: `density` — Risk-Neutral Density (2D Lines)
**Endpoint:** `GET /api/csv/density`
**Nav accent:** `#6bcb77`

**Chart:** Plotly `scatter` traces (2D, NOT 3D)
Three lines: `k_rbs` (blue `#4d96ff`, width 2.5), `k_gauss` (green `#6bcb77`, dash `dash`), `k_cauchy` (red `#ff6b6b`, dash `dot`)

Add:
- `shape: 'tozeroy'` fill for `k_rbs` (light blue fill, opacity 0.15)
- Two vertical `shapes` (dashed yellow `#ffd93d`) at the min and max `x` values where `inside_cone === 0` — these are the cone boundaries
- One `shape` rectangle (red, opacity 0.06) spanning the outside-cone x regions — the "forbidden zone"

Layout: dark bg, x-axis label "Log-return x = ln(S_T/K)", y-axis label "Density", legend top-right.

Include a **toggle button** above the chart: "Linear / Log scale" — switches `yaxis.type` between `'linear'` and `'log'` via `Plotly.relayout`.

---

### Panel 3: `greeks` — Delta and Gamma Surfaces (3D × 2)
**Endpoint:** `GET /api/grid/greeks`
**Nav accent:** `#ffd93d`

**Chart:** Two Plotly `surface` charts side by side (CSS grid `grid-cols-2`)

Left — Delta:
```javascript
{ type: 'surface', x: data.x, y: data.y, z: data.surfaces.delta,
  colorscale: 'Plasma', colorbar: { title: 'Δ' } }
```
Scene: xaxis "Strike K", yaxis "τ_rel" (log), zaxis "Delta", zrange [0, 1]

Right — Gamma:
```javascript
{ type: 'surface', x: data.x, y: data.y, z: data.surfaces.gamma,
  colorscale: 'Viridis', colorbar: { title: 'Γ' } }
```
Scene: xaxis "Strike K", yaxis "τ_rel" (log), zaxis "Gamma"

Both: `camera: { eye: { x: 1.8, y: -1.8, z: 1.0 } }`, dark layout.

---

### Panel 4: `american` — American Put Premium (3D Surface)
**Endpoint:** `GET /api/grid/american`
**Nav accent:** `#ff6b6b`

**Chart:** Single Plotly `surface`
```javascript
{ type: 'surface', x: data.x, y: data.y, z: data.surfaces.premium_pct,
  colorscale: 'Hot', colorbar: { title: 'Premium (%)' } }
```
Scene: xaxis "Strike K", yaxis "τ_rel" (log), zaxis "Early-Exercise Premium (%)"
Camera: `{ eye: { x: -1.5, y: 1.5, z: 1.0 } }`

---

### Panel 5: `barrier` — Barrier Discount (3D Surface)
**Endpoint:** `GET /api/grid/barrier`
**Nav accent:** `#c77dff`

**Chart:** Single Plotly `surface`
```javascript
{ type: 'surface', x: data.x, y: data.y, z: data.surfaces.discount_pct,
  colorscale: 'RdPu', colorbar: { title: 'Discount (%)' } }
```
Scene: xaxis "Barrier Level H", yaxis "τ_rel" (log), zaxis "Discount from Vanilla (%)"
Camera: `{ eye: { x: 1.5, y: -1.5, z: 1.2 } }`

---

### Panel 6: `calibration` — RMSE Landscape (2D Line)
**Endpoint:** `GET /api/csv/calibration`
**Nav accent:** `#ff9a3c`

**Chart:** Plotly `scatter` (2D, NOT 3D)

Traces:
1. RMSE curve: `{ x: data.map(r=>r.tau_rel), y: data.map(r=>r.rmse_bp), type:'scatter', mode:'lines', line:{color:'#4d96ff', width:2} }`
2. Green fill below 1 bp: filter rows where `rmse_bp < 1`, use `fill:'tozeroy'`, `fillcolor:'rgba(107,203,119,0.15)'`
3. Reference line: vertical dashed yellow at `x = 0.05` (τ_true)

Layout:
- `xaxis: { title: 'τ_rel', type: 'log' }` — **Plotly handles log X natively** (no manual mapping needed)
- `yaxis: { title: 'RMSE (basis points)' }`
- Annotation text at x=0.05: "τ_true = 0.05"

---

### Panel 7: `comparison` — Market vs BS vs RBS (3D Triple Surface) ⭐
**Endpoint:** `GET /api/grid/comparison`
**Nav accent:** `#00b4d8`

**Chart:** Three overlaid Plotly `surface` traces on the same axes — **the centrepiece plot**

```javascript
[
  {
    type: 'surface',
    name: 'Market (parametric skew)',
    x: data.x, y: data.y, z: data.surfaces.market_iv,
    colorscale: [[0,'#ff3333'],[1,'#ff9999']],
    opacity: 0.9,
    showscale: false,
    hovertemplate: 'Market IV: %{z:.2f}%<extra></extra>'
  },
  {
    type: 'surface',
    name: 'Black-Scholes (flat)',
    x: data.x, y: data.y, z: data.surfaces.bs_iv,
    colorscale: [[0,'#3333ff'],[1,'#9999ff']],
    opacity: 0.5,
    showscale: false,
    hovertemplate: 'BS IV: %{z:.2f}%<extra></extra>'
  },
  {
    type: 'surface',
    name: 'RBS (frown)',
    x: data.x, y: data.y, z: data.surfaces.rbs_iv,
    colorscale: [[0,'#33aa33'],[1,'#99ff99']],
    opacity: 0.75,
    showscale: false,
    hovertemplate: 'RBS IV: %{z:.2f}%<extra></extra>'
  }
]
```

Layout:
```javascript
{
  scene: {
    xaxis: { title: 'Strike K' },
    yaxis: { title: 'τ_rel', type: 'log' },
    zaxis: { title: 'Implied Volatility (%)', range: [16, 25] },
    camera: { eye: { x: 1.8, y: -1.8, z: 0.9 } }
  },
  legend: { x: 0.02, y: 0.98, bgcolor: 'rgba(0,0,0,0.5)' },
  paper_bgcolor: '#030712', font: { color: '#d1d5db' }, margin: { t: 30 }
}
```

Add **three legend pills** above the chart (not Plotly legend — styled HTML):
- 🔴 Market — realistic skew (OTM puts expensive)
- 🔵 Black-Scholes — flat surface (no smile)
- 🟢 RBS — frown (ATM > wings)

---

## Shared Components

### `<AnnotationCard>` props: `{ title, body, key_insight }`
```jsx
<div className="mt-6">
  <h2 className="text-xl font-semibold text-white mb-3">{title}</h2>
  <div className="inline-flex items-center gap-2 px-3 py-1.5 rounded-full
                  bg-amber-900/30 border border-amber-700/50 text-amber-300
                  text-sm mb-4">
    ⚡ {key_insight}
  </div>
  <div className="bg-gray-900 rounded-xl p-6 border border-gray-800">
    <p className="text-xs text-gray-500 uppercase tracking-wider mb-3">Analysis</p>
    {body.split('\n\n').map((para, i) => (
      <p key={i} className="text-gray-300 text-sm leading-relaxed mb-3">{para}</p>
    ))}
  </div>
</div>
```

### `<PlotlyChart>` wrapper
```jsx
// On mount: Plotly.newPlot(ref.current, traces, layout, { responsive: true, displayModeBar: true })
// On data change: Plotly.react(...)
// Config: { displayModeBar: true, modeBarButtonsToRemove: ['sendDataToCloud'] }
// All 3D charts: add { scrollZoom: false } to prevent accidental zoom
```

### Loading state
```jsx
<div className="flex flex-col items-center justify-center h-96 gap-3">
  <div className="w-8 h-8 border-2 border-blue-500 border-t-transparent rounded-full animate-spin"/>
  <p className="text-gray-400 text-sm">Computing with rel_bs.exe...</p>
</div>
```

### Error state
```jsx
<div className="rounded-lg bg-red-950/50 border border-red-800 p-4 text-red-300 text-sm">
  Failed to load data. Is the server running?
  <code className="block mt-1 text-xs text-red-400">uvicorn server:app --port 8000</code>
</div>
```

---

## Sidebar nav items
```
1. IV Surface      #4d96ff   📈
2. Density         #6bcb77   🔔
3. Greeks          #ffd93d   Δ
4. American Puts   #ff6b6b   ⏱
5. Barrier         #c77dff   🚧
6. Calibration     #ff9a3c   🎯
7. Comparison ⭐   #00b4d8   ⚖
```
Panel 7 has a ⭐ star indicator — it is the centrepiece.

---

## Model parameters card (sidebar bottom)
Fetched from `GET /api/meta` on load:
```
┌─────────────────────┐
│ Model Parameters    │
│ S  = 100  (spot)    │
│ K  = 100  (strike)  │
│ T  = 1 yr           │
│ r  = 5%             │
│ σ  = 20%            │
└─────────────────────┘
```
Small text, gray, monospace values.

---

## Tech stack
- Vite + React + TypeScript
- Tailwind CSS
- **plotly.js-dist** (full build — needed for 3D surface support)
- No state management library — `useState` + `useEffect` per panel with caching in a context
- `fetch` API, no axios
