# Lesson 5 — Wrinkles: Paper 3 and the Markstein correction

Prereqs: Lessons 1–4. Reference: Paper 3 (Hong, Shinar, Fedkiw 2007).

## The problem Paper 3 fixes

Papers 1–2 give a **too smooth** flame front. Real flames are never smooth —
they wrinkle, form cusps, and develop cellular "honeycomb" patterns. Two
physical instabilities are responsible:

### 1. Darrieus–Landau (DL) instability — purely hydrodynamic

From the density jump alone, a planar flame is **linearly unstable**. Any
tiny perturbation grows. The linear growth rate is

```
σ_DL(k) = k · S · [ √( (ρ_f/ρ_h)(1 + ρ_f/ρ_h − ρ_h/ρ_f) ) − 1 ] / (1 + ρ_f/ρ_h)
```

for perturbation wavenumber `k`. Larger `ρ_f/ρ_h` → more unstable. Growth is
proportional to `k`, so *every scale* the grid can represent should grow
at once; smaller scales just grow faster.

In the Paper 1/2 solver this instability is **numerically damped** by the
level-set reinitialization and HJ-WENO advection — both have built-in
dissipation. Result: artificially smooth fronts.

### 2. Thermal-diffusive instability — cells

When the Lewis number `Le = thermal diffusivity / mass diffusivity < 1`
(hydrogen–air, lean methane–air), a bulge in the front pointing into the
fuel gets:
- **Hotter faster** (mass diffuses in more easily than heat escapes),
- **Fuel-richer** (fuel flows in faster than heat can consume it),

so the bulge burns faster and grows. You get cellular patterns — a
honeycomb of convex and concave cells on the front. (For `Le > 1` the front
is stable against cells.)

## The fix: curvature-dependent flame speed

The standard turbulent combustion result (Markstein 1951) is that a curved
flame front burns at a modified speed:

```
S(κ) = S₀ · (1 − M · κ)
```

where:
- `κ = ∇·n = Δφ` is the mean curvature of the front (signed; positive when
  the front bulges into the fuel).
- `M` is the **Markstein length** — a material property with units of length.

### Sign conventions, carefully

- `κ > 0`: front convex toward fuel (a bulge poking into unburned gas).
- `M > 0`: bulges burn **slower** → flat front is stable, perturbations
  decay. Physically: Le > 1, heat escapes too fast for the bulge to keep up.
- `M < 0`: bulges burn **faster** → flat front is unstable, perturbations
  amplify. Physically: Le < 1, mass in faster than heat out.

For a nice wrinkled look you want `M < 0` but not so negative that small
perturbations blow up faster than the nonlinear saturation (from DL
geometry) can contain them. The paper picks `M` on the order of a few grid
cells.

## What has to change in your code

Literally every place `S` appears, substitute `S(κ)`:

### Step 3 — level-set advection

```python
# was:
w = V_f - S * n
# becomes:
kappa = laplacian(phi)              # since |∇φ|=1, κ = Δφ
S_eff = S0 * (1 - M * kappa)
w = V_f - S_eff * n
```

### Step 8 — pressure projection

The jump conditions from Lesson 2 had `S` in them. Now they're curvature-
dependent:

```
[V·n]  = (ρ_f/ρ_h − 1) · S(κ)
[p]    = ρ_f S(κ)² (1 − ρ_f/ρ_h)
```

When you assemble the GFM Poisson RHS, evaluate `S(κ)` at each
front-crossing face using the local `κ` there.

That is the **entirety** of the Paper 3 modification, in its simplest form.
One line of physics, propagated through three places in the code.

## Two-scale amplification (optional)

Even with `M < 0`, a grid-resolved DL instability grows slowly and starts
from numerical noise — you might need thousands of steps to see wrinkles.
Paper 3's contribution beyond `S(κ)` is a **two-scale model** where a
small-scale perturbation `η` on the front evolves under a
Michelson–Sivashinsky-style PDE capturing DL growth + nonlinear saturation,
and the main solver sees `φ_eff = φ − η`.

This amplifies wrinkling on demand without needing infinite grid
resolution. Skip on a first implementation — you'll get serviceable
wrinkles just from `S(κ)`.

## Numerical gotchas

- `κ = Δφ` is very noisy near cusps and away from the front. **Only
  compute and apply it in a narrow band** around `|φ| < c·h`.
- Clamp `κ` to avoid `S(κ)` going negative or unreasonably large. A common
  trick: `S_eff = S0 · max(1 − M·κ, 0.1)`.
- If you use HJ-WENO for advection, it already diffuses short wavelengths;
  you may want to switch to WENO-5 or semi-Lagrangian BFECC for the
  level-set step to give DL room to grow.

## Verification

Set up a 2-D simulation:
- Square domain.
- Initialize `φ` as a horizontal line (half fuel below, half hot above).
- Add a tiny `10⁻³` sinusoidal perturbation to the initial line.
- Run with `M = 0` → perturbation stays flat (numerical damping).
- Run with `M = −0.5·h` → perturbation grows, front develops cusps.

If that works, you've implemented Paper 3.

## Exercise

Starting from the 9-step loop in Lesson 4, annotate each step with "changed
by Paper 3 / unchanged". Convince yourself that steps 1, 2, 4–7, 9 are
untouched; only 3 and 8 see `S → S(κ)`.
