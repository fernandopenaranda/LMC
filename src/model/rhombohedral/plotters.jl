

function LMC.abcplotbandsk(kpoints, p; ylims = [-1, 1], res = (500, 550))
    fig = Figure(size = res)
    LMC.abcplotbands!(fig, real.(abck_bands(kpoints, p)), kpoints, ylimits = ylims, color = :black)
    fig
end

function LMC.abcplotbandsk(fig::Figure, kpoints, p; ylims = [-1, 1], res = (500, 550))
    LMC.abcplotbands!(fig, real.(abck_bands(kpoints, p)), kpoints, ylimits = ylims, color = :orange)
    fig
end

function LMC.abcplotbandsk(ax::Axis, kpoints, p; ylims = [-1, 1], res = (500, 550), color = :black)
    LMC.abcplotbands!(ax, real.(abck_bands(kpoints, p)), ylimits = ylims, color = color)
end


function LMC.abcNplotbandsk(N, kpoints, p; ylims = [-1, 1], res = (500, 550))
    fig = Figure(size = res)
    LMC.abcplotbands!(fig, real.(abck_Nbands(N, kpoints, p)), kpoints, ylimits = ylims, color = :black)
    fig
end

function LMC.abcNplotbandsk(N, kpoints, p, n; ylims = [-1, 1], res = (500, 550))
    fig = Figure(size = res)
    LMC.abcplotbands!(fig, n, real.(abck_Nbands(N, kpoints, p)), kpoints, ylimits = ylims, color = :black)
    fig
end

function LMC.abcNplotbandsk(fig::Figure, N, kpoints, p; ylims = [-1, 1], res = (500, 550))
    LMC.abcplotbands!(fig, real.(abck_Nbands(N, kpoints, p)), kpoints, ylimits = ylims, color = :orange)
    fig
end

function LMC.abcNplotbandsk(ax::Axis, N, kpoints, p; ylims = [-1, 1], res = (500, 550), color = :black)
    LMC.abcplotbands!(ax, real.(abck_Nbands(N, kpoints, p)), ylimits = ylims, color = color)
end

LMC.abcplotbandsk(mat::Array, nk; kw...) = LMC.abcplotbands!(Figure(size = (500, 550)), real.(mat), nk; kw...)

function LMC.abcplotbands!(f::Figure, mat, nk; dots = false, color = missing, ylimits = missing, xlimits = missing)
    ax = Axis(f[1, 1]; xlabel = "k", ylabel = "E [meV]", xlabelsize= 27, ylabelsize= 27, xticklabelsize = 27, yticklabelsize = 27)
    xarr = collect(1:size(mat,2))
    num_points = size(mat,2)
        for i in 1:size(mat, 1)
            lines!(ax, xarr , (mat[i,:]), color = ifelse(isa(color,Missing), :lightgray, color))
        end
 
    ylims!(ax, -0.04*390,0.04*390)
    ax.xticks = ([nk], ["K_ξ"])
    return f
end


function LMC.abcplotbands!(f::Figure, n, mat, nk; dots = false, color = missing, ylimits = missing, xlimits = missing)
    ax = Axis(f[1, 1]; xlabel = "k", ylabel = "E [meV]", title = "ν = $(sum(n))", xlabelsize= 27, ylabelsize= 27, xticklabelsize = 27, yticklabelsize = 27)
    xarr = collect(1:size(mat,2))
    num_points = size(mat,2)
        for i in 1:size(mat, 1)
            lines!(ax, xarr , (mat[i,:]), color = ifelse(isa(color,Missing), :lightgray, color))
        end
 
    ylims!(ax, -0.04*390,0.04*390)
    ax.xticks = ([nk], ["K_ξ"])
    return f
end

function LMC.abcplotbands!(ax::Axis, mat; dots = false, color = :black, ylimits = missing, xlimits = missing)
    xarr = collect(1:size(mat,2))
    num_points = size(mat,2)
    count = 0
        for i in 1:size(mat, 1)
            lines!(ax, xarr , (mat[i,:]), color = color)
            count +=1
        end
    ylims!(ax, -0.05*390,0.05*390)
    ax.xticks = ([div(length(xarr),2)], ["K_ξ"])
end
#_______________________________________________________________________________________
""" computes and plot bands givern a vector of Params_rhombohedral presets """
function LMC.spinfull_plotbands(N, ps; points = 100)
    fig = Figure()
    ax = Axis(fig[1:2, 0]; xlabel = "kx", ylabel = "E [meV]")
    LMC.spinfull_plotbands!(ax, N, ps, points = points)
    return fig
end


function LMC.half_metal_dos(N, p, Delta_Ez, Valley_asym, μlist; evals = 100, η = 0.05)
    ps = half_metal_presets(p, Delta_Ez, Valley_asym)
    ω, js = spinfull_dos(N, ps, η, evals)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "μ [meV]", ylabel = "DOS (a.u.)")
    LMC.plot_dos!(ax, ω, js)
    fig
end

function LMC.quarter_metal_dos(N, p, Delta_Ez, Valley_asym, μlist; evals = 100, η = 0.05)
    ps = quarter_metal_presets(p, Delta_Ez, Valley_asym)
    ω, js = spinfull_dos(N, ps, η, evals)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "μ [meV]", ylabel = "DOS (a.u.)")
    LMC.plot_dos!(ax, ω, js)
    fig
end

function LMC.plot_dos!(ax, ω, js)
    colors = [:black, :black, :gray, :gray]
    styles = [:solid, :dash, :solid, :dash]
    labels = ["+↑", "+↓", "-↑", "-↓"]
    for (i,j) in enumerate(js)
        lines!(ax, ω, j, color = colors[i], linestyle = styles[i], label = labels[i])
    end
    axislegend(ax)
end

function LMC.plotmap(kx, ky, Zs; colrange = missing, colmap = missing, 
    Ω_contr = true, omm_contr = true, fermi_surface = false, with_shift = true, points = 10)
    fig = Figure(size=(1.23*600,600))
    if fermi_surface == false
        ax = Axis(fig[1, 1], xlabel="kx", ylabel="ks", title = "Contributions: Ω: $(Ω_contr), OMM: $(omm_contr), δσijk: $(with_shift)")
    else
        ax = Axis(fig[1, 1], xlabel="kx", ylabel="ks", title = "Fermi surface")
    end
    LMC.plotmap!(ax, kx, ky, Zs)
    cb = Colorbar(fig[1, 2], hm)
end

function LMC.plotmap!(ax, kx, ky, Zs; colmap = missing, colrange = missing)
    if isa(colmap, Missing)
        cmap = cgrad([:black,:red,:white]) 
    else
        cmap = cgrad(colmap) 
    end
    if isa(colrange, Missing)
        hm = heatmap!(ax, kx, ky, real(Zs), rasterize = true)#, colormap = cgrad([:red,:black]))
    else
        hm = heatmap!(ax, kx, ky, real(Zs), colormap = cmap, colorrange= colrange, rasterize = true)
    end
    return hm
    # cb = Colorbar(fig[1, 2], hm)
end

function LMC.plot_characters(N, U, J, Ezlist, νlist, ns)
    characters = character(ns)
    xvect = νlist
    fig = Figure()
    ax= Axis(fig[1,1], ylabel = "D [meV]", xlabel = "ν", title = "Character: N = $(N), U = $(U), J = $(J)")

    scatter!(ax, [100], [100], color = :gray, label = "Symmetric", marker = :rect, markersize = 20)
    scatter!(ax, [100], [100], color = :blue, label = "Half Metal SP", marker = :rect, markersize = 20)
    scatter!(ax, [100], [100], color = :purple, label = "Half Metal VP", marker = :rect, markersize = 20)
    scatter!(ax, [100], [100], color = :lightblue, label = "Quarter Metal/Insulator", marker = :rect, markersize = 20)
    xlims!(ax, xvect[1], xvect[end])
    ylims!(ax, Ezlist[1], Ezlist[end])
    axislegend(ax, position = :rt)
    levels = [-1, 0, 1, 2]
    colors = [:purple, :gray, :blue, :lightblue]
    hm =  heatmap!(ax,  xvect, Ezlist, characters', colormap = colors, levels = levels, colorrange = [-1,2])#, rasterize = true)
    return fig
end