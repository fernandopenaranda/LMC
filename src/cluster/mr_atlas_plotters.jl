include("/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Projects/HeavyFermion_Optics/src/mr_read_atlas")
function plot_mc_atlas(PID::String, B = 10, ϵ = 1e3; ylims = (-5,5))
    ns, σxxx, σxx, rat = reshapemrdata(PID)
    # l = readpresets(PID)
    perm = sortperm(ns)
    nsperm = ns[perm]
    σxxxperm = σxxx[perm]
    σxxperm = σxx[perm]
    ratperm = rat[perm]
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "ν", ylabel = "σ_xx [e^2/h]")#, title = PID * ": GS: " *l.gs[1] * " , evals = $(l.evals[1]), T = $(l.T[1])K")
    lines!(ax, [n for n in nsperm], [ϵ*B*s for s in σxxxperm], label = "σxxx * "*string(B)*"T"*"*"*string(ϵ), color = :red)
    lines!(ax, [n for n in nsperm], [s for s in σxxperm], label = "σxx", color = :black)
    axislegend(ax, position = :lt) 
    ylims!(ax, ylims)
    fig
end

function twoplot_mc_atlas(PID::String, PID2::String, B = 10, ϵ = 1e3; ylims = (-5,5))
    ns, σxxx, σxx, rat = reshapemrdata(PID)
    ns2, σxxx2, σxx2, rat = reshapemrdata(PID2)    
    perm = sortperm(ns)
    perm2 = sortperm(ns2)
    nsperm = ns[perm]
    nsperm2 = ns2[perm2]
    σxxxperm = σxxx[perm]
    σxxxperm2 = σxxx2[perm2]
    σxxperm = σxx[perm]
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "ν", ylabel = "σ_xx [e^2/h]")#, title = PID * ": GS: " *l.gs[1] * " , evals = $(l.evals[1]), T = $(l.T[1])K")
    lines!(ax, [n for n in nsperm], [ϵ*B*s for s in σxxxperm], label = "VH: σxxx * "*string(B)*"T"*"*"*string(ϵ), color = :blue)
    lines!(ax, [n for n in nsperm2], [ϵ*B*s for s in σxxxperm2], label = "QAH: σxxx * "*string(B)*"T"*"*"*string(ϵ), color = :red)
    lines!(ax, [n for n in nsperm], [s for s in σxxperm], label = "σxx", color = :gray)
    axislegend(ax, position = :rb) 
    ylims!(ax, ylims)
    fig
end


function plot_mc_atlas_rat(PID::String, B = 10, ϵ = 1e2; ylims = (-5,5))
    ns, σxxx, σxx, rat = reshapemrdata(PID)
    # l = readpresets(PID)
    
    perm = sortperm(ns)
    nsperm = ns[perm]
    σxxxperm = σxxx[perm]
    σxxperm = σxx[perm]
    ratperm = rat[perm]

    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "ν", ylabel = "σxxx * ["*string(B)*"T]/σ_xx [%]")
    lines!(ax, [n for n in nsperm], [ϵ*B*s for s in σxxxperm ./ (σxxperm+σxxxperm)], color = :black)
    # lines!(ax, [n for n in nsperm], [ϵ*B*s for s in rat], color = :greenq )

    # axislegend(ax, position = :lt) 
    ylims!(ax, ylims)
    fig

end


function plot_in_butterfly(PID::String, n)
    ns, σxxx, σxx, rat = reshapemrdata(PID)
    # l = readpresets(PID)
    perm = sortperm(ns)
    nsperm = ns[perm]
    σxxxperm = σxxx[perm]
    σxxperm = σxx[perm]
    ratperm = rat[perm]
    pos = findclosest(n, nsperm)
    println(nsperm[pos])
    return  σxxxperm[pos]
end

findclosest(val, list) = argmin(abs.(list .- val)) 

E(t) = sin(2π * t)
B(t) = 10*sin(2π*2*t)

function butterfly(ν, PID)
    tlist = collect(0:0.01:1)
    c = plot_in_butterfly(PID, ν)
    return tlist, [B(t)*E(t)*c for t in tlist]
end 

function butterflyvst(ν::Number, PID; color = :black)
    t,val =butterfly(ν, PID)
    fig = Figure()
    # ax = Axis(fig[1,1],ylabel = "j", xlabel = "t")
    # lines!(ax, [i for i in t], [v for v in val])
    ax = Axis(fig[1,1],ylabel = "j(t)", xlabel = "E(t)")
    lines!(ax, [E(i) for i in t], [v for v in val], label = "ν = $(ν)", color = color)

    return fig, ax
end

function butterflyvst!(ax, ν::Number, PID; color = :black)
    t,val =butterfly(ν, PID)
    # lines!(ax, [i for i in t], [v for v in val])
    lines!(ax, [E(i) for i in t], [v for v in val], label = "ν = $(ν)", color = color)
    fig
end

function butterflyvst(ν::Array, PID; ylims = (-0.001,0.001))
    colors = cgrad(:Spectral, length(ν), categorical = true)
    fig, ax = butterflyvst(ν[1], PID, color = colors[1])
    [butterflyvst!(ax, ν[i], PID, color = colors[i]) for i in 2:length(ν)]
    ylims!(ax, ylims)
    Colorbar(fig[1, 2], limits = (minimum(ν), maximum(ν)), size = 10, colormap = colors, vertical = true, ticks=ν, label = "ν")
    fig
end
