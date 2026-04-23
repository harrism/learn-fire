# learn-fire

Working notes and interactive lessons for reimplementing the Fedkiw-group
fire simulation.

## Reading order

1. `lessons/00_overview.md` — the three papers as one story.
2. `lessons/01_level_sets.ipynb` — signed-distance fields, `|∇φ|=1`, reinit. **Interactive.**
3. `lessons/01b_hj_weno5.ipynb` — better advection: HJ-WENO5 + TVD-RK3. **Interactive.**
4. `lessons/02_jump_conditions.md` — the three jumps across the flame front.
5. `lessons/03_ghost_fluid_method.ipynb` — 1D GFM Poisson, including variable-ρ. **Interactive.**
6. `lessons/04_fire_loop.md` — the full Paper 2 per-step algorithm.
7. `lessons/05_markstein_wrinkles.md` — Paper 3 curvature correction.
8. `lessons/06_self_test.md` — questions + answers.

Papers live in `papers/` (download script in `scripts/fetch_papers.sh`).

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
jupyter lab
```

Then open any `.ipynb` under `lessons/`.
