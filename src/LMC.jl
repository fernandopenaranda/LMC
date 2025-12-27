module LMC
    using Arpack
    using LinearAlgebra
    using Cubature
    using ProgressMeter
    using Base.Threads
    using Distributed
    using Dierckx
    using PhysicalConstants
    using PhysicalConstants.CODATA2018
    using Unitful
    using SparseArrays
    using StaticArrays
    using Statistics
    using Roots
    using Parameters
    using QuadGK
    using Optim
    using NLsolve, Interpolations
    using JLD2, CSV, DataFrames
    using Optics_in_the_length_gauge #mine # ] add https://github.com/fernandopenaranda/Optics_in_the_length_gauge

    const kB = (PhysicalConstants.CODATA2018.k_B |> u"eV/K").val
    const μB = (PhysicalConstants.CODATA2018.BohrMagneton |> u"eV/T").val
    const ħ = PhysicalConstants.CODATA2018.ħ
    const e = PhysicalConstants.CODATA2018.e
    const C = ((e^3 / ħ^2) |> u"μA/V^2/s").val
    const C_cd = ((e^2/ħ) |> u"μA/V").val
    const ħ_ev_s = (ħ |> u"eV*s").val
    const ang_to_m = 1e-10

    proj_folder = normpath(joinpath(@__DIR__, "..")) * "src"
    rhomb_folder = proj_folder * "/model/rhombohedral/"
    common_folder = proj_folder * "/comfunctions/"
    figures_folder = proj_folder * "/Paper/Figures/"

    include(proj_folder * "/plot_init.jl")
    include(rhomb_folder * "model.jl")
    include(rhomb_folder * "observables.jl")
    include(rhomb_folder * "filling.jl")
    include(rhomb_folder * "wrapper_lmc.jl")
    include(rhomb_folder * "comp_structs.jl")
    include(rhomb_folder * "spontaneous_sym_breakingmodels.jl")
    include(rhomb_folder * "local_hartree_optimization.jl")
    include(common_folder * "separate_contributions.jl")
    include(proj_folder * "/cluster/cluster_tools.jl")

    export Params_rhombohedral, params_rhombohedral, xxx_lmc_presets, 
        lmc_presets, lmcshift_presets, lmcnoshift_presets, 
        xx_drude_presets, σxyahe_presets, qah_presets
    export Interpolated_params, Computation_params
    export abc_Nlayer, dhxNlg, dhyNlg, dhxxNlg, rzNlg, abc_pentalayer
    export abck_Nbands, abck_bands, abcbz_path_gamma_k_m_gamma
    export rh_filling
    export c_dos, kresolvedlmc, k_Omegain, k_Omegaz, k_d_OMM
    export half_metal_presets, quarter_metal_presets
    export half_metal_plotbands, quarter_metal_plotbands, spinfull_plotbands, spinfull_plotbands!
    export half_metal_dos, quarter_metal_dos, spinfull_dos, plot_dos!, plot_characters
    export μ_α, μ_αT
    export character, reshape_densities, interpolated_dos_ns_Ez, polarization
    export Emin_nαs, Emin_μαs, opt_μs, interpolated_dos, interpolated_n
    export klmc, kresolved_Ωz, kresolved_Ωin, kresolved_dOMM, evalmat, plotmap, plotmap!
    export save_to_csv
    export script_path, slurm_submit_interpolations
end