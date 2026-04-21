# Lesson 0 — Overview

## The three papers as one story

| # | Citation | Role |
|---|----------|------|
| 1 | Nguyen, Fedkiw, Kang 2001 (JCP 172) | **Numerical foundation.** How to run an incompressible Navier–Stokes solver with a sharp discontinuity (the flame front) inside it. |
| 2 | Nguyen, Fedkiw, Jensen 2002 (SIGGRAPH) | **Graphics application.** Uses (1) to build a full fire simulator + blackbody renderer. This is the one you will reimplement end-to-end. |
| 3 | Hong, Shinar, Fedkiw 2007 (SIGGRAPH) | **Realism upgrade.** The front in (1)/(2) is too smooth; this adds Darrieus–Landau and thermal-diffusive instabilities (wrinkles, cells). |

## Mental model

Fire is a **thin reaction zone** (the blue flame sheet) separating two
incompressible gases:

- **Fuel** — cold, dense, ahead of the front. Gets consumed.
- **Hot products** — light, buoyant, behind the front. Glows, carries soot.

The flame front is modeled as a **2-D surface** moving through a 3-D
incompressible flow. Fuel is converted to product at a **flame speed** `S`
(normal to the front, in the fuel frame). Because hot products are ~8× less
dense than fuel, they **expand**: mass conservation across the front forces
a jump in the normal velocity of `(ρ_f/ρ_h − 1)·S`. That expansion is the
physical driver of fire moving outward and upward.

## Fields you'll simulate

| Field | Meaning | Lives where |
|-------|---------|-------------|
| `u` | velocity | MAC grid (staggered), faces |
| `φ` | signed distance to flame front | cell centers |
| `p` | pressure | cell centers |
| `T` | temperature | cell centers |
| `ρ_s` | soot density | cell centers |

Sign convention: `φ < 0` in fuel, `φ > 0` in hot products, `φ = 0` on the
flame front.

## What's actually new vs. a smoke solver

If you already have a Stam/Fedkiw-style smoke simulator, the fire simulator
adds only three modules:

1. **Level-set advection + reinitialization** (Lesson 1).
2. **Jump conditions** relating `u`, `p`, `ρ` across the front (Lesson 2).
3. **GFM modifications to the pressure Poisson solve** (Lesson 3).

Everything else (buoyancy, vorticity confinement, semi-Lagrangian advection,
pressure projection structure) comes straight from the 2001 smoke paper.

## Where we're headed

By the end of these lessons you should be able to:

- Write a level-set module from scratch (Lesson 1 has working code).
- Recite and apply the three jump conditions (Lesson 2).
- Modify a standard CG/ICPCG pressure solver to use GFM (Lesson 3, with a
  1-D working demo).
- Sketch the full per-step loop and know exactly which steps are inherited
  from smoke (Lesson 4).
- Understand Paper 3's one-line modification that produces wrinkles
  (Lesson 5).
