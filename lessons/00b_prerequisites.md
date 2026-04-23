# Lesson 0b — Prerequisites: upwinding and CFL

If you're comfortable with both, skip to Lesson 1. These two ideas sit
underneath every advection step in the fire solver; the fire lessons
assume them as given.

---

## 1. Upwinding

**Rule**: when discretizing an advection term `u·∇φ`, reach toward the
neighbors the information is coming **from**, not where it's going.

### Why

Advection transports values along *characteristics*. For the 1-D equation

```
∂φ/∂t + u φ_x = 0        (u constant for now)
```

solutions are constant along lines of slope `1/u` in `(x, t)`. If `u > 0`,
every value at `(x_i, t+Δt)` came from `(x_i − u·Δt, t)` — to its **left**.

```
        t+Δt    ⋅    ⋅    ⋅    ⋅   φ_i^{n+1}
                       ╲                  ↑
                        ╲                 │  "comes from"
                         ╲                │
        t       ⋅    ⋅    ∙    ⋅    ⋅     │
                     x_{i-1}  x_i         (x_i − u·Δt, t)
                          └── upwind ──┘
```

With `u > 0`, `x_{i-1}` is the upwind neighbor; with `u < 0`, `x_{i+1}` is.

### The three 1st-order finite differences

| Name | Formula | Uses |
|------|---------|------|
| Backward | `(φ_i − φ_{i-1})/h` | `i-1, i` |
| Forward  | `(φ_{i+1} − φ_i)/h` | `i, i+1` |
| Central  | `(φ_{i+1} − φ_{i-1})/(2h)` | `i-1, i+1` |

Von Neumann stability analysis of `∂φ/∂t + u φ_x = 0` with forward Euler:

- `u > 0`, **backward** → stable under CFL (§2). ✓
- `u > 0`, **forward** → **unconditionally unstable**. ✗
- **Central** → **unconditionally unstable** with forward Euler. ✗

The unstable choices reach *toward* where the flow is going, trying to
build the future from data that hasn't physically arrived. Small errors
grow geometrically.

### The stable pattern, in one line of NumPy

```python
phi_x = np.where(u > 0, (phi_i - phi_left)/h, (phi_right - phi_i)/h)
```

This `np.where(u > 0, ...)` pattern appears in every advection function
in the fire lessons.

### Why not just always use central differences — they're 2nd order?

Central differences are fine for **diffusion** (`φ_xx`): information
spreads symmetrically, there's no direction to preserve. For **advection**,
information has a direction. A central stencil has no mechanism to damp
errors that travel *against* the flow, so explicit Euler + central =
unconditionally unstable. (Implicit time integration can rescue it, but
then you're solving a linear system every step.)

### How the idea generalizes

- **Higher-order upwind** (HJ-WENO5, Lesson 1b): reach further in the
  upwind direction. Still upwind-biased; the "5th order" comes from a
  smart weighted combination of candidate stencils.
- **Flux upwinding** (Godunov, Roe): for conservation laws, the flux
  at face `i+½` is evaluated from the upwind state.
- **Semi-Lagrangian**: the CFL-escape cousin. Instead of picking a
  finite-difference stencil, trace the characteristic backward by
  `u·Δt` and interpolate. Geometric upwinding.

### Where it shows up in the fire solver

- `φ` advection: HJ-WENO5 with upwind selection per dimension.
- `T`, `ρ_s` advection: semi-Lagrangian (geometric upwinding).
- `u*` advection: semi-Lagrangian or MacCormack, both upwind-respecting.
- Pressure projection (`Δp = RHS`): elliptic — no direction, central
  differences are the *right* choice.

Rule of thumb: **upwind transport, center diffusion.** Fire needs both.

---

## 2. CFL (Courant–Friedrichs–Lewy) condition

**Rule**: your numerical stencil must "see" at least as far in one step
as the physical wave can travel.

### The intuition

An advection equation carries information at speed `u`. In one timestep
`Δt` a signal travels `u·Δt`. A first-order upwind stencil reaches one
cell `h`. If `u·Δt > h`, the physical signal races past the stencil in a
single step, and the scheme is trying to predict the future from data
that couldn't have reached it yet. Any such scheme goes unstable.

```
physical info cone:   ╲   ╱      numerical stencil:  ░░█░░
                       ╲ ╱                            ↑ reaches one cell
                        ∙
                  (grid point)

CFL satisfied:  physical cone ⊆ stencil cone    → stable
CFL violated:   physical cone ⊃ stencil cone    → unstable (any scheme)
```

This is a **necessary** condition for stability (and for simple explicit
schemes, also sufficient).

### The formula

For 1-D advection with first-order upwind + forward Euler, the **Courant
number** is

```
       u · Δt
C  =  ────────   ≤ 1
          h
```

In multiple dimensions:

```
Δt ≤ h / max(|u|, |v|, |w|)      (grid-aligned CFL)
Δt ≤ h / max(|u| + |v| + |w|)    (stricter; use when safe is required)
```

Advanced schemes trade max C for accuracy:

| Scheme                        | safe max C |
|-------------------------------|-----------:|
| Upwind-1 + forward Euler      |       1.0  |
| HJ-WENO5 + TVD-RK3            |      ~0.6 (use 0.4) |
| Semi-Lagrangian               |   ∞ (unconditionally stable) |

### In practice

You rarely know `max|u|` ahead of time. Each step, compute it and rescale:

```python
u_max = max(np.abs(u).max(), np.abs(v).max(), np.abs(w).max())
dt = C_safety * h / (u_max + 1e-12)   # C_safety ≈ 0.4 for WENO5+RK3
```

For fire, buoyancy and flame expansion make `u` grow as the sim heats up,
so `Δt` shrinks over time. Adaptive timestepping is standard.

### Why CFL is hard

- **Not about accuracy** — it's stability. Violating CFL by 1% doesn't
  give 1% more error; values go to `inf`/`NaN` in a handful of steps.
- **Total compute cost** scales badly. Refining `h → h/2` in 3-D:
  `8×` cells, `2×` more steps (CFL), **16×** total work.
- **Stiffer PDE terms** tighten the bound:

  | Term            | Restriction       | Scaling |
  |-----------------|-------------------|---------|
  | Advection `u·∇` | `Δt ≤ h / |u|`    | `O(h)`  |
  | Diffusion `ν Δ` | `Δt ≤ h² / (2ν)` | `O(h²)` ← often the bottleneck |

### Escapes

- **Implicit time integration** (backward Euler, Crank–Nicolson):
  CFL-free, but each step solves a linear system. Used for diffusion
  terms routinely.
- **Semi-Lagrangian advection**: stable for any `Δt`, but big steps lose
  energy and smear sharp features. Used for `T` and `ρ_s` in Paper 2.
- **IMEX splitting**: implicit on the stiff term, explicit on the rest.

### Where it shows up in the fire solver

- **Level-set advection (Lesson 1b)**: `Δt = 0.4 · h / max|V_f − S·n|`.
  Note the `S·n` term — the front moves at the *effective* velocity, not
  just the fluid velocity.
- **Semi-Lagrangian T, ρ_s transport**: no CFL in principle, but in
  practice you still bound the macro `Δt` by the advection CFL so the
  whole simulation stays consistent.
- **Pressure projection**: elliptic, no CFL.
- **Buoyancy / vorticity confinement forces**: no direct CFL, but they
  raise `|u|` which tightens the advection CFL.

Rule of thumb: **one CFL per explicit term, take the smallest, safety
factor ≈ 0.4.**
