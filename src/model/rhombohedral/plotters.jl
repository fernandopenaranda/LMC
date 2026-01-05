

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
    fig = Figure(size = (600,600))
    ax= Axis(fig[1,1], ylabel = "D [meV]", xlabel = "ν", title = "Phase diagram: N = $(N), U = $(U), J = $(J)")

    scatter!(ax, [100], [100], color = :gray, label = "Symmetric", marker = :rect, markersize = 20)
    scatter!(ax, [100], [100], color = :blue, label = "SPHM/SPHI", marker = :rect, markersize = 20)
    scatter!(ax, [100], [100], color = :purple, label = "VPHM/SPHI", marker = :rect, markersize = 20)
    scatter!(ax, [100], [100], color = :lightblue, label = "QM/QI", marker = :rect, markersize = 20)
    xlims!(ax, xvect[1], xvect[end])
    ylims!(ax, Ezlist[1], Ezlist[end])
    axislegend(ax, position = :ct, orientation = :horizontal, labelsize = 13,
    titlesize = 13,
    padding = (4, 4, 4, 4),
    patchsize = (12, 12))
    levels = [-1, 0, 1, 2]
    colors = [:purple, :gray, :blue, :lightblue]
    hm =  heatmap!(ax,  xvect, Ezlist, characters', colormap = colors, levels = levels, colorrange = [-1,2], rasterize = true)#, rasterize = true)
    return fig
end

function LMC.plot_phasediagrams(pdpath, pdpresetpath)
    ppd = load(pdpresetpath)
    pddata = load(pdpath)
    vals, vals2 = ppd["presets"]
    _, pdpresets = vals
    _, first_value = first(pddata["merged"])
    keylist = keys(first_value)
    Ezs = []
    nss = []
    for (path, d) in pddata["merged"]
        Ez = d["Ezs"][1]
        ns = d["nss"][1]
        push!(Ezs, Ez)
        push!(nss, ns)
    end
    sorted_inds = sortperm(Ezs)
    νlist = first_value["nu_list"]
    plot_characters(pdpresets.N, pdpresets.U, pdpresets.J, Ezs[sorted_inds], νlist, nss[sorted_inds])
end

function LMC.plot_cluster_bandsanddos(pdpath, pdpresetpath, Ezval, νval)
    ppd = load(pdpresetpath)
    pddata = load(pdpath)
    vals, vals2 = ppd["presets"]
    _, pdpresets = vals
    _, first_value = first(pddata["merged"])
    keylist = keys(first_value)
    Ezs = []
    muss = []
    count = 1
    for (path, d) in pddata["merged"]
        Ez = d["Ezs"][1]
        mus = d["mus"][1]
        push!(Ezs, Ez)
        push!(muss, mus)
    end
    νlist = first_value["nu_list"]
    idxnu = argmin(abs.(νlist .- νval))
    idxEz = argmin(abs.(Ezs .- Ezval))
    spinfull_plotbandsanddos(pdpresets.N, pdpresets.p, νlist[idxnu], muss[idxEz][idxnu], Ezs[idxEz], collect(-10:0.1:10); points = 100, evals = 5e3, η = 0.1)
end

function LMC.spinfull_plotbandsanddos(N, p, nαs, μαs, Ez, μlist; points = 100, evals = 100, η = 0.05)

    ps =[Params_rhombohedral(p, μ = μαs[1], ξ = 1, Delta_Ez = Ez, Valley_asym = 0),
    Params_rhombohedral(p, μ = μαs[2], ξ = -1, Delta_Ez = Ez, Valley_asym = 0),
    Params_rhombohedral(p, μ = μαs[3], ξ = 1, Delta_Ez = Ez, Valley_asym = 0),
    Params_rhombohedral(p, μ = μαs[4], ξ = -1, Delta_Ez = Ez, Valley_asym = 0)]

    
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "kx", ylabel = "E [meV]", 
        title = "N = $(N), ν = $(round(sum(nαs),digits = 2)), Ez = $(round(Ez,digits = 2)) meV")
    # spinfull_plotbands!(ax, N, p, μαs, Ez, points = points)
    
    
    LMC.abcNplotbandsk(ax, N, points, ps[1]; ylims = [-1, 1])#, style = :solid)
    LMC.abcNplotbandsk(ax, N, points, ps[2]; ylims = [-1, 1])#, style = :dash)
    LMC.abcNplotbandsk(ax, N, points, ps[3]; ylims = [-1, 1])#, style = :solid)
    LMC.abcNplotbandsk(ax, N, points, ps[4]; ylims = [-1, 1])#, style = :dash)
    
    ax2 = Axis(fig[1,2], xlabel = "μ [meV]", ylabel = "DOS (a.u.)", title = "ν = $(round(sum(nαs),digits = 2)), Ez = $(round(Ez,digits = 2)) meV")
    ω, js = spinfull_dos(N, ps, η, evals, μlist)
    plot_dos!(ax2, ω, js)
    fig
end


function spinfull_plotbands(N, p, nαs, μαs, Ez; points = 100)
    fig = Figure()
    ax = Axis(fig[1:2, 0]; xlabel = "kx", ylabel = "E [meV]", title = "N = $(N), ν = $(round(sum(nαs),digits = 2)), Ez = $(round(Ez,digits = 2)) meV")
    spinfull_plotbands!(ax, N, p, μαs, Ez, points = points)
    return fig
end
"""
Now 1 → K↑, 2 → K'↑, 3 → K↓, 4 → K'↓.
"""
function spinfull_plotbands!(ax, N, p, μαs, Ez; points = 100)
    ps =[Params_rhombohedral(p, μ = μαs[1], ξ = 1, Delta_Ez = Ez, Valley_asym = 0),
        Params_rhombohedral(p, μ = μαs[2], ξ = -1, Delta_Ez = Ez, Valley_asym = 0),
        Params_rhombohedral(p, μ = μαs[3], ξ = 1, Delta_Ez = Ez, Valley_asym = 0),
        Params_rhombohedral(p, μ = μαs[4], ξ = -1, Delta_Ez = Ez, Valley_asym = 0)]

    LMC.abcNplotbandsk(ax, N, points, ps[1]; ylims = [-1, 1])#, style = :solid)
    LMC.abcNplotbandsk(ax, N, points, ps[2]; ylims = [-1, 1])#, style = :dash)
    LMC.abcNplotbandsk(ax, N, points, ps[3]; ylims = [-1, 1])#, style = :solid)
    LMC.abcNplotbandsk(ax, N, points, ps[4]; ylims = [-1, 1])#, style = :dash)
end

function spinfull_dos(N, ps, η, evals, μlist)
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
    labels = ["K↑", "K'↑", "K↓","K'↓"]
    for (i,j) in enumerate(js)
        lines!(ax, ω, j, color = colors[i], linestyle = styles[i], label = labels[i])
    end
    axislegend(ax)
end

LMC.plot_drude(pdpath, pdpresetpath; kws...) = 
    plot_obs(pdpath, pdpresetpath, LMC.plot_drude; kws...)


LMC.plot_ahe(pdpath, pdpresetpath; kws...) = 
    plot_obs(pdpath, pdpresetpath, LMC.plot_ahe; kws...)

LMC.plot_lmc(pdpath, pdpresetpath; kws...) = 
    plot_obs(pdpath, pdpresetpath, LMC.plot_lmc; kws...)

LMC.plot_lmcspin(pdpath, pdpresetpath; kws...) = 
    plot_obs(pdpath, pdpresetpath, LMC.plot_lmcspin; kws...)



function plot_obs(pdpath, pdpresetpath, func; kws...)
    ppd = load(pdpresetpath)
    pddata = load(pdpath)
    vals = ppd["presets"]
    interpPID = vals["comp_struct"].interpPID
    pdPID = vals["comp_struct"].phasediagPID
    PID = vals["comp_struct"].PID
    evals = vals["comp_struct"].evals
    _, first_value = first(pddata["merged"])
    keylist = keys(first_value)
    Ezs = []
    obsmat = []
    for (path, d) in pddata["merged"]
        Ez = d["Ezs"][1]
        o = d["obs_list"]
        push!(Ezs, Ez)
        push!(obsmat, o)
    end
    sorted_inds = sortperm(Ezs)
    νlist = first_value["nu_list"]
    mat = reshape_observables(obsmat[sorted_inds])
    if func == LMC.plot_ahe || func == LMC.plot_lmc
        valley_ordering = reorder_valleys(pdpath, pdpresetpath)
        matt = abs.(mat) .* valley_ordering
    elseif func == LMC.plot_drude || func == LMC.plot_lmcspin
                valley_ordering = ones(size(mat,1), size(mat,2))
        matt = mat
    end
    func(PID, pdPID, interpPID, evals, Ezs[sorted_inds], νlist, matt; kws...)
end


#_______________________________________________________________________________________
# DEALING WITH DEGENERACIES
#_______________________________________________________________________________________
# after the inclusion of Hunds there is only Valley degeneracy
# which sets the sign of the qah with filling 
# we arbitrarily choose the following filling sequence for the valleys
""" 
the spin degeneracy is lifted by the hunds coupling that makes, however
valley degeneracy is still present. One can choose among the different filling
orderings in the valley. This will result in different hall conductivities"""
function reorder_valleys(pdpath, pdpresetpath)
    ppd = load(pdpresetpath)
    pddata = load(pdpath)
    vals = ppd["presets"]
    Ezs = []
    nss = []
    for (path, d) in pddata["merged"]
        Ez = d["Ezs"][1]
        ns = d["nss"][1]
        push!(Ezs, Ez)
        push!(nss, ns)
       
    end
    sorted_inds = sortperm(Ezs)
    sorted_nss = nss[sorted_inds]
    dim1 = size(Ezs,1)
    dim2 = size(nss[1],1)
   
    qah_sign_mat = zeros(Float64, dim1, dim2)
    for i in 1:dim1
        for j in 1:dim2-1
            qah_sign_mat[i,j] = v_filling_order(sorted_nss[i][j])
        end
    end
    return qah_sign_mat
end

""" this order will corresponde to filling first the K valley and
then the K' of arbitrary spin if negative or and the opposite for possitive
fillings, this is completely arbitrary, the same sign could be observed
everywhere """
function v_filling_order(ns)
    ifelse(sum(ns) >0, 1, -1)
end


function reshape_observables(obs)
    dim1 = length(obs)
    dim2 = length(obs[1])
    mat = zeros(Float64, dim1, dim2)
    for i in 1:dim1
        for j in 1:dim2-1
            mat[i,j] = obs[i][j]
        end 
    end
    return mat
end

#_______________________________________________________________________________________
LMC.plot_lmc(PID::String, pdPID, interpPID, evals, Ezlist, νlist, mat; kws...) = 
    plot_mat(PID, pdPID, interpPID, evals, Ezlist, νlist, mat, "LMC", "σxxx [1/T*e²/h]"; kws...)

LMC.plot_lmc(N::Number, U, J, Ezlist, νlist, mat; kws...) = 
    plot_mat(N, U, J, Ezlist, νlist, mat, "LMC", "σxxx [1/T*e²/h]"; kws...)

LMC.plot_lmcspin(PID::String, pdPID, interpPID, evals, Ezlist, νlist, mat; kws...) = 
    plot_mat(PID, pdPID, interpPID, evals, Ezlist, νlist, mat, "LMC spin", "σxxx [1/T*e²/h]"; kws...)

LMC.plot_lmcspin(N::Number, U, J, Ezlist, νlist, mat; kws...) = 
    plot_mat(N, U, J, Ezlist, νlist, mat, "LMC spin", "σxxx [1/T*e²/h]"; kws...)
#_________________________________________
#_______________________________________________________________________________________
# Drude
#_________________________________________________________________________________________
LMC.plot_drude(PID::String, pdPID, interpPID, evals, Ezlist, νlist, mat; kws...) = 
    plot_mat(PID, pdPID, interpPID, evals, Ezlist, νlist, mat, "σxx", "σxx [e²/h]"; kws...)

LMC.plot_drude(N::Number, U, J, Ezlist, νlist, mat; kws...) = 
    plot_mat(N, U, J, Ezlist, νlist, mat, "σxx", "σxx [e²/h]"; kws...)
#_________________________________________________________________________________________
# AHE. xy
#_________________________________________________________________________________________
LMC.plot_ahe(PID::String, pdPID, interpPID, evals, Ezlist, νlist, mat; kws...) = 
    plot_mat(PID, pdPID, interpPID, evals, Ezlist, νlist, mat, "AHE", "σxy [e²/h]"; kws...)

LMC.plot_ahe(N::Number, U, J, Ezlist, νlist, mat; kws...) = 
    plot_mat(N, U, J, Ezlist, νlist, mat, "AHE", "σxy [e²/h]"; kws...)

#_________________________________________________________________________________________
# Plotter
#_________________________________________________________________________________________
function plot_mat(N::Number, U, J, Ezlist, νlist, matt, label, labellegend; colorrange = [-maximum(abs.(matt)), maximum(abs.(matt))])
    fig = Figure()
    ax= Axis(fig[1,1], xlabel = "Ez [meV]", ylabel = "ν", title = "$(label): N = $(N), U = $(U), J = $(J)")
    ylims!(ax,  νlist[1], νlist[end])
    xlims!(ax, Ezlist[1], Ezlist[end])
    colors = cgrad([:lightblue,:black,:red])
    # hm =  heatmap!(ax, Ezlist,  νlist, mat, colormap = colors, colorrange = [-maximum(mat), maximum(mat)])#, colorrange = [-1,2])#, rasterize = true)
    hm =  heatmap!(ax, νlist, Ezlist, matt', colormap = colors,  rasterize = true, colorrange = [-1,2])#, rasterize = true)
    cb = Colorbar(fig[1,2], hm, label = labellegend)
    return fig
end

function plot_mat(PID::String, pdPID, interpPID, evals, Ezlist, νlist, matt, label, labellegend; colorrange = [-maximum(abs.(matt)), maximum(abs.(matt))])
    fig = Figure(size=(600,600))
    ax= Axis(fig[1,1], xlabel = "ν", ylabel = "Ez [meV]", title = "$(label): PID = $(PID), intPID = $(interpPID), pdPID = $(pdPID), evals = $(evals)")
    xlims!(ax,  νlist[1], νlist[end])
    ylims!(ax, Ezlist[1], Ezlist[end])
    colors = cgrad([:lightblue,:black,:red])
    # hm =  heatmap!(ax, Ezlist,  νlist, mat, colormap = colors, colorrange = [-maximum(mat), maximum(mat)])#, colorrange = [-1,2])#, rasterize = true)
    hm =  heatmap!(ax, νlist, Ezlist, matt', colormap = colors, colorrange = colorrange, rasterize = true)#, colorrange = [-1,2])#, rasterize = true)
    cb = Colorbar(fig[1,2], hm, label = labellegend)
    return fig
end