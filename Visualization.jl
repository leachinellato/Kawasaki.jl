module Visualization

using Plots

export animate_positions, plot_observables

"""
    animate_positions(frames, kT, J; fps=30, name="positions.gif")

Animate a sequence of spin configurations.
"""
function animate_positions(frames::AbstractVector, kT::Real, J::Real;
                           fps::Int = 30, name::String = "positions.gif")
    my_palette = palette([:yellow, :blue])
    anim = @animate for i in eachindex(frames)
        heatmap(frames[i];
                c = my_palette,
                clims = (-1, 1),
                colorbar_title = "Spin",
                aspect_ratio = :equal,
                title = "J = $J, kT = $kT, step = $i")
    end
    gif(anim, name; fps = fps)
end

"""
    plot_observables(result; filename=nothing)

Plot energy and magnetization traces from `run_simulation`.
"""
function plot_observables(result; filename::Union{Nothing,String} = nothing)
    p1 = plot(result.energy;
              xlabel = "snapshot",
              ylabel = "E",
              legend = false,
              title = "Energy")
    p2 = plot(result.magnetization;
              xlabel = "snapshot",
              ylabel = "M",
              legend = false,
              title = "Magnetization (should be constant)")
    p = plot(p1, p2; layout = (2, 1), size = (700, 500))
    filename !== nothing && savefig(p, filename)
    return p
end

end # module
