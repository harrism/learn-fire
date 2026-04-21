# Learning the Fedkiw Fire Papers

## Context

You are preparing to reimplement a physically-based fire simulation. The three
papers below form a tight lineage from Ron Fedkiw's group at Stanford:

1. **Nguyen, Fedkiw, Kang 2001 (JCP)** — the *numerical foundation*: how to
   advance an incompressible flow that has a sharp flame discontinuity in it.
2. **Nguyen, Fedkiw, Jensen 2002 (SIGGRAPH)** — the *graphics application*:
   uses (1) to build a full fire simulator plus a blackbody renderer.
3. **Hong, Shinar, Fedkiw 2007 (SIGGRAPH)** — the *realism upgrade*: the
   flame-front model in (1)/(2) is too smooth; this adds the Darrieus–Landau
   and thermal-diffusive instabilities that produce wrinkles and cellular
   structure on real flames.

Read them in this order. (2) is the one you will reimplement end-to-end; (1)
supplies the key subroutines; (3) is an optional "v2" layer once (2) runs.

This document teaches the physics and numerics you need to understand each
paper. A follow-up plan will turn this into an implementation plan once we
agree on language, scope, and rendering ambitions.

---

## Shared mental model

Fire is a *thin reaction zone* (the blue flame sheet) separating two
incompressible gases:

- **Cold fuel / premixed gas** — ahead of the front. Gets consumed.
- **Hot gaseous products** — behind the front. Buoyant, glowing, carries soot.

The reaction zone is a few hundred microns thick in reality, effectively a
2-D surface at any graphics scale. Fuel is converted to product at a
characteristic *flame speed* `S` (normal to the front, measured in the fuel
frame). Because hot products are much less dense, they *expand*: mass is
conserved across the front, so the normal velocity must jump. That jump is
the hard part — a naive incompressible solver smears it out.

Core data structures (same across all three papers):

- `u` — velocity on a MAC grid (staggered), incompressible (`∇·u = 0` on each
  side of the front, with a specific jump across it).
- `φ` — level set, signed distance to the flame front. `φ < 0` fuel, `φ > 0`
  hot products.
- `T` — temperature (advected scalar, cooled by radiation).
- `Y` / `ρ_soot` — density of soot / smoke product (advected scalar).
- `p` — pressure (solved via Poisson each step).

Constants you will hard-code:
- `ρ_f` fuel density, `ρ_h` hot-product density (`ρ_f/ρ_h` ≈ 8 is typical).
- `S` flame speed (normal speed of front relative to fuel, e.g. 0.5 m/s).
- `T_ignition`, `T_max`, `T_air`; cooling coefficient `c_T`.
- Buoyancy `α` (temperature) and `β` (smoke weight).

---

## Paper 1 — Nguyen, Fedkiw, Kang 2001 (JCP 172)

### The problem it solves

"How do I run an incompressible Navier–Stokes solver when the fluid has a
discontinuity in density and normal velocity across a moving interface?"

Before this paper, graphics-style fluid solvers (Stam's Stable Fluids) assume
one continuous fluid. A flame front violates that: ρ jumps, u·n jumps, p
jumps. You need to *capture* those jumps inside the pressure projection
rather than smear them.

### The jump conditions (memorize these)

Let `n` be the unit normal from fuel into hot product, `D` the front's
geometric speed, `V_f` and `V_h` the fluid velocities on each side, `S` the
flame speed (rate at which fuel is consumed, measured relative to the fuel).

Mass conservation across the moving front:

```
ρ_f (V_f - D) · n = ρ_h (V_h - D) · n = ρ_f · S   (= mass flux through front)
```

From this you derive:

```
V_f · n = D · n + S                      (fuel flows into front at speed S)
V_h · n = D · n + (ρ_f / ρ_h) · S
[V · n] = V_h·n − V_f·n = (ρ_f/ρ_h − 1) · S      ← the key expansion jump
```

So the *hot side* moves faster outward than the *cold side* — by a constant
multiplier of `S`. That is the physical origin of fire "pushing" outward.

Front advection (your level set update):

```
∂φ/∂t + w · ∇φ = 0,   w · n = D · n = V_f · n − S     (advect from the fuel side)
```

Pressure jump (from momentum balance, viscous terms dropped):

```
[p] = (ρ_f − ρ_h)(D · n)^2 ≈ ρ_f (1 − ρ_f/ρ_h) · S^2   (small; often 0 in 2002)
```

### The Ghost Fluid Method (GFM)

Instead of discretizing across the jump (which smears), at each grid cell
near the interface you:

1. Keep the *real* fluid's value on one side.
2. Build a *ghost* value on the other side that makes the stencil continuous
   by adding the known jump `[·]`.
3. Use standard finite differences on each side, with ghosts filling in.

For the pressure Poisson solve, this gives a symmetric positive-definite
system whose RHS encodes `[p]` and `[u·n]`. The resulting `∇p` correction,
when subtracted from the intermediate velocity, enforces `∇·u = 0` on each
side *and* the correct velocity jump across the front. That's the whole
trick.

### What you take from this paper

- Jump conditions above (write them on a sticky note).
- Level-set advection with velocity `V_f − S·n` (or the hot-side equivalent).
- Modified pressure Poisson with GFM for the flame front.
- Extension-velocity / fast-sweeping reinitialization of `φ` to keep it a
  signed distance function.

---

## Paper 2 — Nguyen, Fedkiw, Jensen 2002 (SIGGRAPH)

### What it adds

This is the fire paper. Paper 1 is the engine; Paper 2 is the car.
Contributions beyond (1):

- Full SIGGRAPH-style simulator: buoyancy, vorticity confinement, temperature
  and soot transport.
- Fuel model: supports both **gaseous premixed fuel** (volume source) and
  **solid fuel** (flame front offset from solid surface by distance `d`).
- Blackbody radiation rendering.

### The per-step algorithm

```
1. Reinitialize φ to signed distance (fast sweeping / PDE-based).
2. Compute flame-front normal n = ∇φ/|∇φ|, curvature κ = ∇·n.
3. Advect φ:  φ^{n+1} = φ^n − Δt (V_f − S·n) · ∇φ   (upwind / HJ-WENO).
4. Advect temperature T and soot density ρ_s:
      semi-Lagrangian or MacCormack, using the correct-side velocity.
5. Cooling:  T ← T − Δt · c_T · ((T − T_air)/T_max)^4     (Stefan–Boltzmann).
6. External forces on hot side:
      f = α (T − T_air) ẑ − β ρ_s ẑ                       (buoyancy).
   Plus vorticity confinement: f_vc = ε h (N × ω),
      ω = ∇×u, N = ∇|ω|/|∇|ω||.
7. Integrate u* = u + Δt (−(u·∇)u + ν∇²u + f).
8. Pressure projection with GFM jump conditions across φ=0:
      solve ∇·(∇p/ρ) = ∇·u*/Δt with jumps [p] and [u·n] from Paper 1.
      u^{n+1} = u* − Δt ∇p/ρ.
9. Rendering pass (offline): ray-march, accumulate blackbody emission
      from temperature, extinction from soot.
```

### Gaseous vs. solid fuel

- **Gaseous**: fuel fills a source region where `φ < 0` and is replenished
  (e.g. a gas jet). Velocity inside the fuel is the incompressible field you
  solve for.
- **Solid**: the fuel doesn't flow. Instead you keep `φ` as a signed distance
  from the flame front, which sits a small offset `d` outside the solid. The
  "fuel side" velocity is prescribed (zero or a puff-off velocity) and only
  the hot side is simulated as a fluid. This is how they do the candles and
  matches.

### Blackbody rendering (essential but often skipped in v1)

Each voxel emits Planck-spectrum radiance `B_λ(T)`. Integrate along camera
rays:

```
L = ∫ exp(−∫ κ ds') · C_bb(T) · ρ_s ds
```

Colors come from `T` via a precomputed `T → RGB` LUT (CIE-weighted Planck).
Soot both emits and absorbs; fuel region is transparent (or faintly blue
from a separate chemiluminescence term they add).

### What you take from this paper

- The full 9-step loop above.
- Buoyancy + vorticity confinement (standard Fedkiw recipe, reused from his
  2001 smoke paper).
- Blackbody LUT for rendering.
- Decision point for your reimplementation: gaseous, solid, or both.

---

## Paper 3 — Hong, Shinar, Fedkiw 2007 (SIGGRAPH)

### Why the previous model looks too smooth

Flames have two intrinsic hydrodynamic instabilities that Paper 2 does not
reproduce:

1. **Darrieus–Landau (DL) instability** — a *purely hydrodynamic* instability
   from the density jump across the front. Any planar flame is linearly
   unstable; perturbations grow and wrinkle it. The simple GFM front is
   stable because numerical diffusion kills the growth.
2. **Thermal-diffusive instability** — when the Lewis number `Le < 1` (fuel
   diffuses faster than heat), bulges toward fuel get hotter faster, burn
   faster, and amplify; you get cellular "honeycomb" patterns on the front.

Real flames show wrinkles, cusps, and cellular cells because of these. The
2002 solver smooths them away.

### The fix: curvature-dependent flame speed (Markstein model)

Replace the constant `S` with:

```
S(κ) = S_0 (1 − M · κ)
```

where `κ` is the front's mean curvature and `M` is the Markstein length.
- `M > 0` stabilizes bulges (convex-toward-fuel parts burn slower).
- `M < 0` amplifies them → cellular instability (what you want for visual
  richness).

Plugging `S(κ)` back into the jump conditions from Paper 1 is the bulk of
the mathematical work: `[u·n]`, `[p]`, and the level-set advection all pick
up `κ`-dependent terms.

### Two-scale front model

Resolving DL on a single uniform grid is expensive, so they:

- Track large-scale front position with a standard level set (as in 2002).
- Layer a *thin-flame model* with a perturbation `η` on the front that
  obeys its own evolution equation capturing DL growth plus nonlinear
  saturation (Michelson–Sivashinsky-like).
- Feed the perturbed front back into the main solver as an effective `φ_eff
  = φ − η`.

This is cheaper than resolving DL volumetrically and controllable via the
Markstein length and a growth-rate parameter.

### Practical extras

- Thermal diffusion term added to the temperature transport (small but
  required for the cellular pattern to appear).
- Careful treatment of cusps — when the front self-intersects, they project
  back to a signed distance with local surgery.

### What you take from this paper

- `S → S(κ)` — a one-line change that already roughens the front.
- Optional: the Michelson–Sivashinsky-style perturbation equation for the
  full wrinkling look.
- A reason *why* a direct reimplementation of Paper 2 looks too waxy, and a
  knob (Markstein length) to tune.

---

## Suggested reading order and reimplementation on-ramp

1. **Stam 1999 "Stable Fluids"** (if rusty) — gives you the semi-Lagrangian
   + pressure-projection pattern everything here builds on.
2. **Fedkiw, Stam, Jensen 2001 "Visual Simulation of Smoke"** — adds
   vorticity confinement and temperature buoyancy. Paper 2 inherits its
   outer loop from this.
3. **Paper 1 (Nguyen–Fedkiw–Kang 2001)** — read §2 (jump conditions) and §3
   (GFM pressure solve) carefully; skim the rest.
4. **Paper 2 (2002)** — this is your spec.
5. **Paper 3 (2007)** — read once, defer implementation until Paper 2
   stands up.

### Open questions before we can write an implementation plan

- Target language/runtime? (Python + NumPy for learning vs. C++/CUDA for
  performance vs. Taichi or WGSL for a middle ground.)
- 2D first (much easier, teaches the ideas) or straight to 3D?
- Fuel type: gaseous jet, solid (candle/match), or both?
- Rendering: matplotlib quicklooks, volumetric ray marcher with blackbody,
  or write out OpenVDB for external rendering?
- Goal for Paper 3: implement now, stub, or skip?

Once we pick these, the follow-up plan file will enumerate modules (grid,
advection, level set, GFM Poisson, combustion, renderer), dependency order,
and per-module verification steps.

## Fetching the PDFs

Plan mode is read-only, so I located the canonical URLs but haven't
downloaded them yet. On approval, I will `mkdir -p papers/` and `curl` the
three PDFs:

- Paper 1 (2001 JCP, flame discontinuities):
  `https://physbam.stanford.edu/~fedkiw/papers/cam2000-19.pdf`
  → `papers/2001_nguyen_fedkiw_kang_flame_discontinuities.pdf`
- Paper 2 (2002 SIGGRAPH, fire):
  `https://physbam.stanford.edu/papers/stanford2002-02.pdf`
  → `papers/2002_nguyen_fedkiw_jensen_fire.pdf`
- Paper 3 (2007 SIGGRAPH, wrinkled flames):
  `https://www.cs.ucr.edu/~shinar/papers/2007_wrinkled_flames.pdf`
  → `papers/2007_hong_shinar_fedkiw_wrinkled_flames.pdf`

## Verification of *this* learning plan

This file is a reference document, not code. To confirm it's useful for you:

- Can you, after reading this, restate the three jump conditions from
  Paper 1 without looking?
- Can you list the 9 steps of the Paper 2 per-frame loop?
- Can you explain in one sentence what the Markstein length controls?

If yes to all three, we're ready to plan the reimplementation.
