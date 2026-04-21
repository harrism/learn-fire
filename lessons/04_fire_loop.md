# Lesson 4 — The full Paper 2 per-step loop

Prereqs: Lessons 1–3. Reference: Paper 2 §§3–6.

## The 9 steps

```
1. Reinitialize φ to signed distance.                            (Lesson 1)
2. Compute n = ∇φ,  κ = Δφ   (valid because |∇φ|=1).              (Lesson 1)
3. Advect φ:
      ∂φ/∂t + (V_f − S n) · ∇φ = 0    (HJ-WENO or semi-Lagrangian)
4. Advect T and ρ_s:
      semi-Lagrangian or MacCormack, using the correct-side velocity.
5. Radiative cooling:
      T ← T − Δt · c_T · ((T − T_air)/T_max)^4
6. External forces on the hot side only:
      f_buoy = α(T − T_air) ẑ − β ρ_s ẑ
      f_vc   = ε h (N × ω),   ω = ∇×u,  N = ∇|ω|/|∇|ω||
7. Integrate:   u* = u + Δt (−(u·∇)u + ν∇²u + f)
8. Pressure projection with GFM jumps across φ=0:                 (Lesson 3)
      solve  ∇·(∇p/ρ) = ∇·u*/Δt   with [p], [V·n] on front-crossing faces
      u^{n+1} = u* − Δt ∇p/ρ
9. Rendering pass (offline): ray-march, blackbody emission + soot extinction.
```

## What's new vs. a smoke solver

| Step | Inherited from smoke (Fedkiw-Stam-Jensen 2001)? | Fire-specific? |
|------|--------------------------------------------------|----------------|
| 1    | —                                                | ✓ level set    |
| 2    | —                                                | ✓ level set    |
| 3    | —                                                | ✓ level set + flame speed offset |
| 4    | ✓ (structure)                                    | side-correct velocity |
| 5    | ✓ verbatim                                       |                |
| 6    | ✓ verbatim                                       |                |
| 7    | ✓ verbatim                                       |                |
| 8    | ✓ (structure)                                    | GFM jumps      |
| 9    | —                                                | ✓ blackbody LUT |

If you have a working smoke solver, implementing fire reduces to:
- Add the level-set module (Lesson 1): ~200 lines.
- Modify the pressure solver to accept jumps (Lesson 3): ~100 lines.
- Add the blackbody renderer: ~100 lines.

Everything else is reused.

## Fuel models

### Gaseous premixed fuel

The fuel fills a source region where `φ < 0` and is replenished each step.
Example: a gas jet nozzle continuously injects fuel at some inflow velocity.
The **fuel side itself is a fluid** — its velocity is part of the simulated
incompressible field.

### Solid fuel (candles, matches)

The fuel doesn't flow. Instead:
- Keep `φ` as the signed distance from the flame *front*, which sits a small
  offset `d` outside the solid surface.
- The "fuel side" velocity near the front is prescribed by a puff-off model
  (some small outward gas velocity from pyrolysis) — not solved.
- Only the hot side is simulated as Navier-Stokes.

Pick one for your reimplementation. Gaseous is conceptually simpler and
produces a good "flamethrower / gas flame" look; solid is what you want for
candles, torches, and burning objects.

## External force details

### Buoyancy

```
f_buoy = α (T − T_air) ẑ − β ρ_s ẑ
```

- First term: hot gas rises. `α ≈ 0.01` per Kelvin (tune visually).
- Second term: soot is heavy, sinks. `β ≈ 0.1` (tune).

Both act **only on the hot side** (zero inside fuel).

### Vorticity confinement

You know this from smoke. Unchanged:

```
ω = ∇ × u
η = ∇ |ω|
N = η / |η|
f_vc = ε · h · (N × ω)       # h is grid spacing; ε ~ 10 for fire
```

`ε` wants to be **larger for fire than smoke** — real flames are more
turbulent than drifting smoke, and confinement is how you recover that
turbulence on coarse grids.

### Radiative cooling

```
T ← T − Δt · c_T · ((T − T_air)/T_max)^4
```

The `^4` is Stefan-Boltzmann. `c_T` sets how fast the plume cools above the
flame — tune so that the tip of your plume has dropped from ~2000 K (yellow)
to ~800 K (dull red) over the visible column height.

## Boundary conditions

- **Solver domain walls**: `u · n = 0` (free-slip) or `u = 0` (no-slip).
  Doesn't matter much if the domain is big enough.
- **Fuel source**: Dirichlet on `u` inside the source region (inflow).
- **Flame front (the interior)**: *this is not a boundary in the usual
  sense* — it's handled by the GFM jumps, not by zeroing anything out.

## Rendering

Each voxel emits Planck-spectrum radiance `B_λ(T)`. Integrate along camera
rays:

```
L = ∫ exp(−∫ κ_ext ds') · C_bb(T) · ρ_s ds
```

Practically:
1. Precompute a `T → RGB` LUT by integrating `B_λ(T) · colormatch(λ)` over
   the visible range for temperatures 0–3000 K. (Use CIE 1931 for physical
   accuracy, or fudge it.)
2. For each ray, march through the grid, accumulate emission and extinction.
3. Fuel side (`φ < 0`) is transparent (or add a faint chemiluminescent blue
   proportional to distance from the front).

Start with matplotlib-style 2D slice output. Upgrade to a proper volumetric
renderer once the physics is right.

## Implementation order (for your reimplementation)

1. Port / write smoke solver first (if you don't have one).
2. Add level-set module. Verify with circle advection test, reinit test.
3. Add GFM pressure solve. Verify with 1-D jump test from Lesson 3, then
   2-D expanding-disc test (fuel disc surrounded by hot products —
   should expand radially at `(ρ_f/ρ_h − 1)·S`).
4. Glue together: full 2-D fire. Visualize `T` with matplotlib — should
   look like a flame silhouette.
5. Add blackbody rendering.
6. Move to 3-D.
7. (Optional) Paper 3 wrinkling — one-line change (Lesson 5).
