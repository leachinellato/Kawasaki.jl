#!/usr/bin/env julia
# Run Kawasaki-Metropolis and produce an animation + observables plot.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "src", "KawasakiMetropolis.jl"))
include(joinpath(@__DIR__, "..", "src", "Visualization.jl"))

using .KawasakiMetropolis
using .Visualization
using Random

# ---------------------------------------------------------------- parameters
const N              = 64
const J              = 1.0
const kT             = 0.5
const n_steps        = 2_000_000
const snapshot_every = 5_000
const m              = 0.0      # magnetization density (conserved)
const seed           = 42

rng = Xoshiro(seed)

# ---------------------------------------------------------------- run
@info "Running Kawasaki-Metropolis" N J kT n_steps

result = @time run_simulation(; N = N, J = J, kT = kT,
                                n_steps = n_steps,
                                snapshot_every = snapshot_every,
                                m = m, rng = rng)

@info "Acceptance ratio" result.acceptance
@info "Final energy per site" result.energy[end] / (N * N)

# ---------------------------------------------------------------- plots
mkpath(joinpath(@__DIR__, "..", "Plots"))
outdir = joinpath(@__DIR__, "..", "Plots")

animate_positions(result.frames, kT, J;
                  fps = 30,
                  name = joinpath(outdir, "animation_kT=$(kT).gif"))

plot_observables(result;
                 filename = joinpath(outdir, "observables_kT=$(kT).png"))

@info "Done. Output written to $outdir"
