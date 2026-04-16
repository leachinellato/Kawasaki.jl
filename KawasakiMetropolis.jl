module KawasakiMetropolis

using Random
using Statistics

export generate_init_cond, generate_init_cond_magnetization,
       total_energy, magnetization,
       kawasaki_step!, metropolis!,
       run_simulation

"""
    generate_init_cond(N; rng=Random.default_rng())

Random ±1 configuration on an N×N lattice.
"""
function generate_init_cond(N::Integer; rng::AbstractRNG = Random.default_rng())
    return rand(rng, Int8[-1, 1], N, N)
end

"""
    generate_init_cond_magnetization(N, m; rng=Random.default_rng())

Configuration with fixed magnetization density `m ∈ [-1, 1]`. Since Kawasaki
dynamics conserves magnetization, this is the usual starting point for
phase-separation studies.
"""
function generate_init_cond_magnetization(N::Integer, m::Real;
                                          rng::AbstractRNG = Random.default_rng())
    @assert -1 <= m <= 1 "magnetization density must be in [-1, 1]"
    total = N * N
    n_up = round(Int, (1 + m) * total / 2)
    spins = fill(Int8(-1), total)
    spins[1:n_up] .= Int8(1)
    shuffle!(rng, spins)
    return reshape(spins, N, N)
end

"""
    total_energy(conf, J)

Full lattice energy with periodic boundary conditions.
H = -J Σ_<ij> s_i s_j
"""
function total_energy(conf::AbstractMatrix{<:Integer}, J::Real)
    N = size(conf, 1)
    E = 0
    @inbounds for j in 1:N, i in 1:N
        s = conf[i, j]
        # Sum only right and down neighbors to avoid double counting
        E += -J * s * (conf[mod1(i + 1, N), j] + conf[i, mod1(j + 1, N)])
    end
    return float(E)
end

"""
    magnetization(conf)

Total magnetization Σ s_i.
"""
magnetization(conf::AbstractMatrix{<:Integer}) = sum(conf)

"""
    local_dE(conf, i, j, di, dj, J)

Energy change for swapping site (i,j) with its neighbor at (i+di, j+dj),
computed from local neighborhoods only. Returns 0 for same-spin pairs.
"""
@inline function local_dE(conf::AbstractMatrix{<:Integer},
                           i::Int, j::Int, di::Int, dj::Int, J::Real)
    N = size(conf, 1)
    s1 = conf[i, j]
    i2, j2 = mod1(i + di, N), mod1(j + dj, N)
    s2 = conf[i2, j2]
    s1 == s2 && return 0.0  # swap is a no-op

    # Neighbors of site 1, excluding site 2
    nsum1 = conf[mod1(i - 1, N), j] + conf[mod1(i + 1, N), j] +
            conf[i, mod1(j - 1, N)] + conf[i, mod1(j + 1, N)] - s2
    # Neighbors of site 2, excluding site 1
    nsum2 = conf[mod1(i2 - 1, N), j2] + conf[mod1(i2 + 1, N), j2] +
            conf[i2, mod1(j2 - 1, N)] + conf[i2, mod1(j2 + 1, N)] - s1

    # ΔE = -J * [(s2 - s1) * nsum1 + (s1 - s2) * nsum2]
    #    = -J * (s2 - s1) * (nsum1 - nsum2)
    return -J * (s2 - s1) * (nsum1 - nsum2)
end

const NEIGHBOR_OFFSETS = ((1, 0), (-1, 0), (0, 1), (0, -1))

"""
    kawasaki_step!(conf, J, kT; rng=Random.default_rng()) -> (accepted, ΔE)

Propose and (possibly) perform a single Kawasaki nearest-neighbor swap in place.
Returns a tuple `(accepted::Bool, ΔE::Float64)` so the caller can track running
energy without recomputation.
"""
function kawasaki_step!(conf::AbstractMatrix{<:Integer}, J::Real, kT::Real;
                        rng::AbstractRNG = Random.default_rng())
    N = size(conf, 1)
    i = rand(rng, 1:N)
    j = rand(rng, 1:N)
    di, dj = NEIGHBOR_OFFSETS[rand(rng, 1:4)]

    ΔE = local_dE(conf, i, j, di, dj, J)

    # Metropolis acceptance
    if ΔE <= 0 || rand(rng) < exp(-ΔE / kT)
        i2, j2 = mod1(i + di, N), mod1(j + dj, N)
        conf[i, j], conf[i2, j2] = conf[i2, j2], conf[i, j]
        return true, ΔE
    end
    return false, 0.0
end

"""
    run_simulation(; N=32, J=1.0, kT=0.5, n_steps=10^6,
                     snapshot_every=1000, m=0.0, rng=Random.default_rng())

Run Kawasaki-Metropolis and return a NamedTuple with:
- `frames`        :: Vector of lattice snapshots (copies, taken every `snapshot_every`)
- `energy`        :: Running energy trace (one entry per snapshot)
- `magnetization` :: Running magnetization trace (should be constant)
- `acceptance`    :: Overall acceptance ratio

The starting configuration is built with magnetization density `m`.
"""
function run_simulation(; N::Integer = 32,
                         J::Real = 1.0,
                         kT::Real = 0.5,
                         n_steps::Integer = 10^6,
                         snapshot_every::Integer = 1000,
                         m::Real = 0.0,
                         rng::AbstractRNG = Random.default_rng())
    conf = generate_init_cond_magnetization(N, m; rng = rng)
    E = total_energy(conf, J)
    M = magnetization(conf)

    n_frames = n_steps ÷ snapshot_every + 1
    frames = Vector{Matrix{Int8}}(undef, 0)
    sizehint!(frames, n_frames)
    energy_trace = Float64[]
    sizehint!(energy_trace, n_frames)
    mag_trace = Int[]
    sizehint!(mag_trace, n_frames)

    push!(frames, copy(conf))
    push!(energy_trace, E)
    push!(mag_trace, M)

    n_accepted = 0
    for step in 1:n_steps
        accepted, ΔE = kawasaki_step!(conf, J, kT; rng = rng)
        if accepted
            E += ΔE
            n_accepted += 1
        end
        if step % snapshot_every == 0
            push!(frames, copy(conf))
            push!(energy_trace, E)
            push!(mag_trace, M)  # conserved by construction
        end
    end

    return (frames = frames,
            energy = energy_trace,
            magnetization = mag_trace,
            acceptance = n_accepted / n_steps)
end

end # module
