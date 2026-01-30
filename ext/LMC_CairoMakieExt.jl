module LMC_CairoMakieExt
    using CairoMakie
    using LMC, JLD2
    println("Loaded plotting extension")
    proj_folder = normpath(joinpath(@__DIR__, "..")) * "src"
    rhomb_folder = proj_folder * "/model/rhombohedral/"
    tbg_folder = proj_folder * "/model/tbg/"
    include(rhomb_folder * "plotters.jl")
    include(tbg_folder * "hf_plotters.jl")
end