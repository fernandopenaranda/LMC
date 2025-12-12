#_________________________________________________________________________________________
# PRESETS FOR A QUARTER AND HALF METAL MODEL
#_________________________________________________________________________________________
"""given a gap size and valley asymmetry params it returns a set of parameters for the computation of a given observable assuming a half metal.
Valley asymmetry is not really required I think but I leave it just in case"""
function half_metal_presets(p, Delta_Ez, Valley_asym = 0)
    pvsu = Params_rhombohedral(p, ξ = 1, Delta_Ez = 0, Valley_asym = Valley_asym)
    pvsd = Params_rhombohedral(p, ξ = 1, Delta_Ez = 0, Valley_asym = Valley_asym)
    nvsu = Params_rhombohedral(p, ξ = -1, Delta_Ez = Delta_Ez, Valley_asym = Valley_asym)
    nvsd = Params_rhombohedral(p, ξ = -1, Delta_Ez = Delta_Ez, Valley_asym = Valley_asym)
    return [pvsu, pvsd, nvsu, nvsd]
end

"""given a gap size and valley asymmetry params it returns a set of parameters for the computation of a given observable assuming a quarter metal."""
function quarter_metal_presets(p, Delta_Ez, Valley_asym = 0)
    pvsu = Params_rhombohedral(p, ξ = 1, Delta_Ez = 0, Valley_asym = Valley_asym)
    pvsd = Params_rhombohedral(p, ξ = 1, Delta_Ez = Delta_Ez, Valley_asym = Valley_asym)
    nvsu = Params_rhombohedral(p, ξ = -1, Delta_Ez = Delta_Ez, Valley_asym = Valley_asym)
    nvsd = Params_rhombohedral(p, ξ = -1, Delta_Ez = Delta_Ez, Valley_asym = Valley_asym)
    return [pvsu, pvsd, nvsu, nvsd]
end

#_________________________________________________________________________________________
# BANDS for quarter and half metals
#_________________________________________________________________________________________
""" bands of the half_metal_model """
function half_metal_plotbands(p::Params_rhombohedral, N, Delta_Ez, Valley_asym = 0; points = 100)
    ps = half_metal_presets(p, Delta_Ez, Valley_asym)
    return spinfull_plotbands(N, ps, points = points)
end

""" bands of the quarter_metal_model """
function quarter_metal_plotbands(p::Params_rhombohedral, N, Delta_Ez, Valley_asym = 0; points = 100)
    ps = quarter_metal_presets(p, Delta_Ez, Valley_asym)
    return spinfull_plotbands(N, ps, points = points)
end

""" computes and plot bands givern a vector of Params_rhombohedral presets """
function spinfull_plotbands(N, ps; points = 100)
    fig = Figure()
    ax = Axis(fig[1:2, 0]; xlabel = "kx", ylabel = "E [meV]")
    spinfull_plotbands!(ax, N, ps, points = points)
    return fig
end

function spinfull_plotbands!(ax, N, ps; points = 100)
    abcNplotbandsk(ax, N, points, ps[1]; ylims = [-1, 1], color = :black)#, style = :solid)
    abcNplotbandsk(ax, N, points, ps[2]; ylims = [-1, 1], color = :black)#, style = :dash)
    abcNplotbandsk(ax, N, points, ps[3]; ylims = [-1, 1], color = :gray)#, style = :solid)
    abcNplotbandsk(ax, N, points, ps[4]; ylims = [-1, 1], color = :gray)#, style = :dash)
end

#_________________________________________________________________________________________
# DOS methods for quarter and half metals
#_________________________________________________________________________________________

function half_metal_dos(N, p, Delta_Ez, Valley_asym, μlist; evals = 100, η = 0.05)
    ps = half_metal_presets(p, Delta_Ez, Valley_asym)
    ω, js = spinfull_dos(N, ps, η, evals)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "μ [meV]", ylabel = "DOS (a.u.)")
    plot_dos!(ax, ω, js)
    fig
end

function quarter_metal_dos(N, p, Delta_Ez, Valley_asym, μlist; evals = 100, η = 0.05)
    ps = quarter_metal_presets(p, Delta_Ez, Valley_asym)
    ω, js = spinfull_dos(N, ps, η, evals)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "μ [meV]", ylabel = "DOS (a.u.)")
    plot_dos!(ax, ω, js)
    fig
end

function spinfull_dos(N, ps, η, evals)
    nps = [xxx_lmc_presets(N, ps[i]) for i in 1:length(ps)]
    ω, j1 = c_dos(nps[1], μlist, η = η, evals = evals)
    ω, j2 = c_dos(nps[2], μlist, η = η, evals = evals)
    ω, j3 = c_dos(nps[3], μlist, η = η, evals = evals)
    ω, j4 = c_dos(nps[4], μlist, η = η, evals = evals)
    return ω, [j1,j2,j3,j4]
end

function plot_dos!(ax, ω, js)
    colors = [:black, :black, :gray, :gray]
    styles = [:solid, :dash, :solid, :dash]
    labels = ["+↑", "+↓", "-↑", "-↓"]
    for (i,j) in enumerate(js)
        lines!(ax, ω, j, color = colors[i], linestyle = styles[i], label = labels[i])
    end
    axislegend(ax)
end
