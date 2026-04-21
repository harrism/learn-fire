# Lesson 6 — Self-test

Try to answer these from memory. Answers below.

## Questions

1. What does `|∇φ| = 1` mean physically, and what three things does it buy
   you cheaply?

2. If `φ` is advected under a velocity field for 50 steps, does its zero
   set still represent the flame front? Does it still satisfy `|∇φ| = 1`?

3. Write the three jump conditions across the flame front. Which are
   spatial constants in the simplest model?

4. Why does mass conservation force a velocity jump of exactly
   `(ρ_f/ρ_h − 1)·S`?

5. If you removed the GFM and used a plain Laplacian for the pressure
   Poisson solve, what would the simulation look like?

6. In the GFM 1-D Poisson problem, where does the jump `[p]` enter: in
   the matrix, in the right-hand side, or both?

7. Which steps in the Paper 2 per-frame loop are inherited verbatim from a
   smoke solver? Which are genuinely new?

8. What is the Markstein length `M` and what sign of `M` produces wrinkles?

9. What single substitution, made in two places, turns a Paper 2
   implementation into a Paper 3 implementation?

10. Why is `κ = Δφ` only true when `|∇φ| = 1`?

---

## Answers

1. **`|∇φ| = 1` means φ is a signed-distance function — moving one unit in
   space changes φ by one unit.** It buys you: (a) unit normal
   `n = ∇φ` for free (no division), (b) mean curvature `κ = Δφ` for free,
   (c) exact distance-to-front `|φ|` at every grid point. It also makes
   the linear interpolation used to locate the front on grid edges
   (`θ = φ_i/(φ_i − φ_{i+1})`) exact.

2. **Zero set: yes.** Advection preserves the zero set regardless.
   **`|∇φ| = 1`: no.** Shear in the velocity field stretches / compresses
   φ, so `|∇φ|` drifts. That's why you reinitialize each step.

3. Sticky note:
   - `[ρ] = ρ_h − ρ_f`  (constant)
   - `[V·n] = (ρ_f/ρ_h − 1)·S`  (constant in simplest model)
   - `[p] = ρ_f S² (1 − ρ_f/ρ_h)`  (constant; often set to 0)

   All three are spatial constants when `ρ_f`, `ρ_h`, `S` are constants.
   Paper 3 turns `S → S(κ)` so they pick up curvature dependence.

4. In the front's frame, mass crosses in from the fuel side at rate `ρ_f·S`
   per unit area. That mass must exit on the hot side at the same rate,
   but the hot side has density `ρ_h`, so it must move at `(ρ_f/ρ_h)·S`
   relative to the front. In the lab frame: `V_f·n = D·n + S`,
   `V_h·n = D·n + (ρ_f/ρ_h)·S`, subtract.

5. The jump would smear over 2–3 cells. The sharp expansion kick that
   drives upward flame motion would weaken. The flame would look like
   buoyant smoke, not fire: puffy, soft, no crisp blue core, no strong
   upward draft from mass-density mismatch.

6. **Right-hand side only.** The LHS Laplacian matrix (with `1/ρ`
   coefficients that differ per side) stays symmetric positive-definite.
   The jumps `[p]` and `[V·n]` appear as extra source terms on the RHS for
   cells adjacent to the front. This is the whole reason GFM is
   attractive — you keep your existing CG/ICPCG solver.

7. **Verbatim from smoke**: advect T and ρ_s (step 4), cooling (5),
   buoyancy (6a), vorticity confinement (6b), u* integration (7),
   projection *structure* (8). **New**: level-set reinit (1), geometry
   computation n, κ (2), level-set advection (3), GFM modifications to
   step 8, blackbody rendering (9).

8. `M` is a material length scale (order of flame thickness) that makes
   the flame speed curvature-dependent: `S(κ) = S₀(1 − M κ)`.
   **`M < 0`** → bulges into fuel burn faster → instability grows → wrinkles.

9. Replace `S` with `S₀(1 − M·κ)`. Two places it appears: (a) level-set
   advection velocity `w = V_f − S·n` in step 3; (b) the jumps `[V·n]`
   and `[p]` used to assemble the GFM Poisson RHS in step 8.

10. Curvature is defined as `κ = ∇·n`, where `n = ∇φ/|∇φ|` is the unit
    normal. Only when `|∇φ| = 1` does `n = ∇φ` directly, making
    `κ = ∇·∇φ = Δφ`. Otherwise you have to compute
    `κ = ∇·(∇φ/|∇φ|)`, which involves extra derivatives and is more
    sensitive to noise.

---

## Ready to implement?

If you got 8/10 on the conceptual answers (details don't matter, structure
does), you're ready. Proceed to:

1. Set up the Python environment (`requirements.txt`).
2. Write a level-set module (SDF + reinit + advection) with unit tests.
3. Write a 1-D GFM Poisson solver and reproduce the notebook in
   Lesson 3.
4. Port or write a 2-D smoke solver if you don't have one.
5. Glue together, run in 2-D, debug visually.
6. Go 3-D once 2-D looks right.

Open a fresh conversation and ask for an implementation plan — we'll
turn these lessons into a concrete module breakdown and dependency order.
