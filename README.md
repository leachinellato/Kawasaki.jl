# Kawasaki-Metropolis on the 2D Ising lattice

Metropolis implementation of Kawasaki (spin-exchange) dynamics for the
2D Ising model on a square lattice with periodic boundary conditions.

Author: Leandro M. Chinellato <chinellato.leandro@gmail.com>

## What it does

Simulates the Ising Hamiltonian

H = -J Σ_<ij> s_i s_j

with local conservation of magnetization: at each step two nearest-neighbor
spins are swapped, accepted with Metropolis probability `min(1, exp(-ΔE/kT))`.
ΔE is computed from the local neighborhoods only (O(1) per step), and the
total energy is tracked incrementally.

## Repository layout

```
.
├── Project.toml           Environment (Julia ≥ 1.9)
├── src/
│   ├── KawasakiMetropolis.jl   Core module: dynamics + observables
│   └── Visualization.jl        Plots/animation helpers
├── scripts/
│   └── run.jl                  Driver: runs a simulation and saves outputs
├── test/
│   └── runtests.jl             Sanity checks
└── Plots/                      Output (git-ignored)
```

## Reproducing

From the repo root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/run.jl
```

This will write `animation_kT=0.5.gif` and `observables_kT=0.5.png` to
`Plots/`. Change `N`, `kT`, `n_steps`, `m` at the top of `scripts/run.jl`.

## Running tests

```bash
julia --project=. test/runtests.jl
```

The tests check:
- ground-state energy of the fully polarized configuration
- magnetization conservation under Kawasaki dynamics
- agreement between the incremental energy and a full recomputation
- that equal-spin swaps have ΔE = 0
- that the local ΔE matches the brute-force full-lattice difference

## What's different from v0.1

- Local ΔE instead of recomputing the full energy every step (orders of
  magnitude speedup).
- In-place swaps; no `deepcopy` in the hot loop.
- Incremental energy tracking, running magnetization, acceptance ratio.
- Configurable initial magnetization (Kawasaki conserves it, so this is the
  physically meaningful control).
- Snapshot every `snapshot_every` steps, with `copy` (not reference) so
  frames are independent.
- Removed unused `MKL` and `LaTeXStrings` imports.
- `Project.toml` committed so the environment is reproducible.
- Code split into a core module and a visualization module; a standalone
  driver script; a test suite.
