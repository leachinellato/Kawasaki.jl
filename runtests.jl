using Test
using Random

include("KawasakiMetropolis.jl")
using .KawasakiMetropolis

@testset "KawasakiMetropolis" begin

    @testset "Initial conditions" begin
        rng = Xoshiro(1)
        c = generate_init_cond(8; rng = rng)
        @test size(c) == (8, 8)
        @test all(x -> x == 1 || x == -1, c)
        @test eltype(c) == Int8

        c2 = generate_init_cond_magnetization(10, 0.4; rng = rng)
        @test sum(c2) == round(Int, 0.4 * 100)
    end

    @testset "Energy consistency (ground state)" begin
        # Fully polarized state on N×N with J=1 has E = -2 N^2
        N = 6
        up = fill(Int8(1), N, N)
        @test total_energy(up, 1.0) ≈ -2 * N^2
        down = fill(Int8(-1), N, N)
        @test total_energy(down, 1.0) ≈ -2 * N^2
    end

    @testset "Magnetization conservation" begin
        rng = Xoshiro(123)
        result = run_simulation(; N = 16, J = 1.0, kT = 1.0,
                                  n_steps = 50_000,
                                  snapshot_every = 1000,
                                  m = 0.25, rng = rng)
        @test all(result.magnetization .== result.magnetization[1])
    end

    @testset "Running energy matches recomputed energy" begin
        # Track energy incrementally and compare to full recomputation
        rng = Xoshiro(7)
        N = 12
        J = 1.0
        kT = 0.8
        conf = generate_init_cond_magnetization(N, 0.0; rng = rng)
        E = total_energy(conf, J)
        for _ in 1:20_000
            accepted, ΔE = kawasaki_step!(conf, J, kT; rng = rng)
            if accepted
                E += ΔE
            end
        end
        @test E ≈ total_energy(conf, J) atol = 1e-8
    end

    @testset "Equal-spin swaps have ΔE = 0" begin
        # Build a config with a definite equal-spin pair
        conf = fill(Int8(1), 5, 5)
        conf[3, 3] = Int8(1)
        conf[3, 4] = Int8(1)
        # Swap (3,3) <-> (3,4); both are +1 so ΔE must be zero
        @test KawasakiMetropolis.local_dE(conf, 3, 3, 0, 1, 1.0) == 0.0
    end

    @testset "ΔE agrees with full energy difference" begin
        # Compare local_dE against brute-force recomputation after a swap
        rng = Xoshiro(99)
        N = 8
        J = 1.3
        conf = generate_init_cond_magnetization(N, 0.0; rng = rng)
        for _ in 1:200
            i = rand(rng, 1:N); j = rand(rng, 1:N)
            di, dj = rand(rng, [(1,0), (-1,0), (0,1), (0,-1)])
            E_before = total_energy(conf, J)
            ΔE = KawasakiMetropolis.local_dE(conf, i, j, di, dj, J)
            i2, j2 = mod1(i+di, N), mod1(j+dj, N)
            conf[i,j], conf[i2,j2] = conf[i2,j2], conf[i,j]
            E_after = total_energy(conf, J)
            @test ΔE ≈ E_after - E_before atol = 1e-10
            # swap back
            conf[i,j], conf[i2,j2] = conf[i2,j2], conf[i,j]
        end
    end

end
