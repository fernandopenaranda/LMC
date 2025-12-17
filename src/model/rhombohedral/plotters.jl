

function abcplotbandsk(kpoints, p; ylims = [-1, 1], res = (500, 550))
    fig = Figure(size = res)
    abcplotbands!(fig, real.(abck_bands(kpoints, p)), kpoints, ylimits = ylims, color = :black)
    fig
end

function abcplotbandsk(fig::Figure, kpoints, p; ylims = [-1, 1], res = (500, 550))
    abcplotbands!(fig, real.(abck_bands(kpoints, p)), kpoints, ylimits = ylims, color = :orange)
    fig
end

function abcplotbandsk(ax::Axis, kpoints, p; ylims = [-1, 1], res = (500, 550), color = :black)
    abcplotbands!(ax, real.(abck_bands(kpoints, p)), ylimits = ylims, color = color)
end


function abcNplotbandsk(N, kpoints, p; ylims = [-1, 1], res = (500, 550))
    fig = Figure(size = res)
    abcplotbands!(fig, real.(abck_Nbands(N, kpoints, p)), kpoints, ylimits = ylims, color = :black)
    fig
end

function abcNplotbandsk(N, kpoints, p, n; ylims = [-1, 1], res = (500, 550))
    fig = Figure(size = res)
    abcplotbands!(fig, n, real.(abck_Nbands(N, kpoints, p)), kpoints, ylimits = ylims, color = :black)
    fig
end

function abcNplotbandsk(fig::Figure, N, kpoints, p; ylims = [-1, 1], res = (500, 550))
    abcplotbands!(fig, real.(abck_Nbands(N, kpoints, p)), kpoints, ylimits = ylims, color = :orange)
    fig
end

function abcNplotbandsk(ax::Axis, N, kpoints, p, n; ylims = [-1, 1], res = (500, 550), color = :black)
    abcplotbands!(ax, n, real.(abck_Nbands(N, kpoints, p)), ylimits = ylims, color = color)
end

abcplotbandsk(mat::Array, nk; kw...) = abcplotbands!(Figure(size = (500, 550)), real.(mat), nk; kw...)

function abcplotbands!(f::Figure, mat, nk; dots = false, color = missing, ylimits = missing, xlimits = missing)
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


function abcplotbands!(f::Figure, n, mat, nk; dots = false, color = missing, ylimits = missing, xlimits = missing)
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

function abcplotbands!(ax::Axis, mat; dots = false, color = :black, ylimits = missing, xlimits = missing)
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