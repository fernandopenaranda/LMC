module LMC_CairoMakieExt
    using CairoMakie
    using LMC
    println("Loaded plotting extension")
    proj_folder = normpath(joinpath(@__DIR__, "..")) * "src"
    rhomb_folder = proj_folder * "/model/rhombohedral/"
    include(rhomb_folder * "plotters.jl")
end