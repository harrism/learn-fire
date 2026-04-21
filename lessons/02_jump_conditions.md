# Lesson 2 — Jump conditions across the flame front

Prereq: Lesson 1 (level sets). Reference: Paper 1 §2.

## Setup

Zoom in on the flame front Γ at one point. Locally it's a plane with unit
normal `n` pointing **from fuel into hot products**. Define:

| Symbol | Meaning |
|--------|---------|
| `D`    | velocity of the geometric front itself (only `D·n` matters) |
| `V_f`  | fluid velocity on the fuel side |
| `V_h`  | fluid velocity on the hot side |
| `ρ_f`  | fuel density (constant) |
| `ρ_h`  | hot-product density (constant, `ρ_h < ρ_f`) |
| `S`    | flame speed — how fast the front eats fuel in the fuel rest frame |

Typical values: `ρ_f/ρ_h ≈ 8`, `S ≈ 0.5 m/s`.

## Jump 1 — normal velocity

### Fuel side

In the fuel's rest frame, the front advances into unburned gas at speed `S`.
In the lab frame: `V_f · n − D · n = S`, so

```
V_f · n = D · n + S                                            (A)
```

Fuel is crossing **through** the front into the hot side at rate `S`.

### Hot side — mass conservation

Mass flux leaving fuel side = mass flux entering hot side (the front
generates no mass, only converts fuel to product):

```
ρ_f (V_f · n − D · n) = ρ_h (V_h · n − D · n)
         ρ_f · S      = ρ_h (V_h · n − D · n)
    V_h · n = D · n + (ρ_f / ρ_h) · S                          (B)
```

### The jump

Subtract (A) from (B):

```
[V · n]  ≡  V_h · n − V_f · n  =  (ρ_f/ρ_h − 1) · S
```

This is **a spatial constant** (for constant `ρ_f`, `ρ_h`, `S`). With
`ρ_f/ρ_h = 8` and `S = 0.5 m/s`, `[V·n] = 3.5 m/s`. That's the expansion
kick — hot gas has to shoot away from the front because mass in must equal
mass out and the outgoing gas is 8× less dense, so it has to move 8× faster.
The *difference* is 7× faster, which is where `(ρ_f/ρ_h − 1)` comes from.

## Jump 2 — pressure

Momentum balance across the front (divergence theorem on a pillbox
straddling Γ, drop viscous terms):

```
[ p + ρ (V · n − D · n)² ] = 0
```

i.e. the quantity inside brackets is continuous across the front. Plug in
`V_f · n − D · n = S` and `V_h · n − D · n = (ρ_f/ρ_h) S`:

```
p_f + ρ_f S²  =  p_h + ρ_h · (ρ_f/ρ_h)² · S²
               =  p_h + (ρ_f² / ρ_h) · S²

[p] ≡ p_h − p_f  =  ρ_f S² − (ρ_f²/ρ_h) S²
                 =  ρ_f S² (1 − ρ_f/ρ_h)
```

This is **negative** (since `ρ_f/ρ_h > 1`): the hot side sits at slightly
lower pressure, which is what sucks fuel across the front. Magnitude is
small — for `ρ_f = 1.0 kg/m³`, `S = 0.5 m/s`, `ρ_f/ρ_h = 8`, you get
`[p] ≈ −1.75 Pa`. The 2002 paper often just sets `[p] = 0` with no visible
difference; we'll follow that simplification.

## Jump 3 — density (trivial)

```
[ρ]  =  ρ_h − ρ_f                (spatial constant)
```

## The sticky note

```
╔═══════════════════════════════════════════════════════╗
║ [ρ]    =  ρ_h − ρ_f                                   ║
║ [V·n]  =  (ρ_f/ρ_h − 1) · S                           ║
║ [p]    =  ρ_f S² (1 − ρ_f/ρ_h)    (often set to 0)    ║
║                                                       ║
║ Front advection:                                      ║
║   ∂φ/∂t + w·∇φ = 0,   w = V_f − S·n                   ║
╚═══════════════════════════════════════════════════════╝
```

All three jumps are **spatial constants** in the simplest model. That is
what makes GFM tractable for fire — you can bake the jumps into the Poisson
RHS once per step, instead of computing them dynamically per face.

## Worked example

A laminar candle flame. Use:
- `ρ_f = 1.0`, `ρ_h = 0.125` (so `ρ_f/ρ_h = 8`).
- `S = 0.5 m/s`.
- The flame front is stationary in the lab (`D · n = 0`) because fuel is
  flowing up into it at speed `S` from the wick.

Then:
- `V_f · n = 0 + 0.5 = 0.5 m/s` (upward, into the front) ✓
- `V_h · n = 0 + 8 · 0.5 = 4.0 m/s` (the hot plume shoots up at 4 m/s) ✓
- `[V·n] = 4.0 − 0.5 = 3.5 m/s` ✓
- `[p] = 1.0 · 0.25 · (1 − 8) = −1.75 Pa`.

Those numbers *are* the simulation, minus buoyancy and vorticity.

## Exercise

The problem: in Lesson 3 we'll implement the GFM Poisson solve. Before you
get there, convince yourself that you could, if asked, evaluate all three
jumps at any point on the front given `ρ_f`, `ρ_h`, `S`. Write them out
without looking back at this page. When you can, move on.
