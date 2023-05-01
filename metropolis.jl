using LinearAlgebra
using MKL
using Random
using Statistics
using Plots
using LaTeXStrings

function animate_positions(frames, kT, J; fps=30, name="positions.gif")
    my_palette = palette([:yellow, :blue])
    anim = @animate for i = 1:length(frames)
        heatmap(frames[i], c = my_palette, colorbar_title = "Spin", title="J = $J, kT = $kT")        
    end
    gif(anim, name, fps=fps)
end

function approve_flip(P)
    if rand() < P
        return "flip"
    else
        return "noflip"
    end
end

function Generate_init_cond(N::Int64)
    #conf = rand([-1,1],N,N)
    conf = rand([Int8(-1), Int16(1)], N,N) 
    return conf
end


function calculate_energy(config::Matrix, J::Float64)
    N = size(config, 1)
    energy = 0
    for i in 1:N
        for j in 1:N
            spin = config[i, j]
            neighbor_sum = config[mod1(i-1, N), j] + config[mod1(i+1, N), j] +
                           config[i, mod1(j-1, N)] + config[i, mod1(j+1, N)]
            energy += -J * spin * neighbor_sum
        end
    end
    return energy / 2 #for avoid double counting
end


function Kawasaki_exachange(conf::Matrix)
    N = size(conf,1)
    ex_indi, ex_indj = rand(1:N), rand(1:N)
    inter = rand(([1,0],[-1,0],[0,1],[0,-1]))
    #println("\n indice a intercambiar: ($ex_indi, $ex_indj)")
    #println("\n itercambia con: $inter")
    new_confg = deepcopy(conf)
    new_confg[ex_indi, ex_indj] = conf[mod1(ex_indi + inter[1],N) , mod1(ex_indj + inter[2],N)]
    new_confg[mod1(ex_indi + inter[1],N) , mod1(ex_indj + inter[2],N)] = conf[ex_indi, ex_indj]
    return new_confg
end

function Metropolis(J::Float64, conf::Matrix, kT::Float64, n_steps::Int64)
    
    frames = []
    for _ = 1:n_steps
        push!(frames, conf)    
        new_conf = Kawasaki_exachange(conf)
        E_i = calculate_energy(conf, J)
        E_f = calculate_energy(new_conf, J)
        ΔE = E_f - E_i
        P = exp(-ΔE/(kT))
        if ΔE <= 0
            conf = new_conf
            new_conf = nothing
        else
            if approve_flip(P) == "flip"
                conf = new_conf
                new_conf = nothing
            else 
                new_conf = nothing
            end
        end
    end
    return frames
end


N = 15
kT = 0.1
init_conf = Generate_init_cond(N)
n_steps = 5000
J = 1.0

frames = Metropolis(J, init_conf, kT, n_steps)



animate_positions(frames, kT,J; fps=90, name="animation_kT="*string(kT)*".gif")