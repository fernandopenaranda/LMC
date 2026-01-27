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
    using Glob
    using Statistics

    const kB = (PhysicalConstants.CODATA2018.k_B |> u"eV/K").val
    const μB = (PhysicalConstants.CODATA2018.BohrMagneton |> u"eV/T").val
    const ħ = PhysicalConstants.CODATA2018.ħ
    const e = PhysicalConstants.CODATA2018.e
    const C = ((e^3 / ħ^2) |> u"μA/V^2/s").val
    const C_cd = ((e^2/ħ) |> u"μA/V").val
    const ħ_ev_s = (ħ |> u"eV*s").val
    const ang_to_m = 1e-10

    const proj_folder = normpath(joinpath(@__DIR__, "..")) * "src"
    const rhomb_folder = proj_folder * "/model/rhombohedral/"
    const tbg_folder = proj_folder * "/model/tbg/"
    const common_folder = proj_folder * "/comfunctions/"
    const figures_folder = proj_folder * "/Paper/Figures/"

    include(proj_folder * "/plot_init.jl")
    include(rhomb_folder * "model.jl")
    include(rhomb_folder * "observables.jl")
    include(rhomb_folder * "filling.jl")
    include(rhomb_folder * "wrapper_lmc.jl")
    include(rhomb_folder * "comp_structs.jl")
    include(rhomb_folder * "spontaneous_sym_breakingmodels.jl")
    include(rhomb_folder * "local_hartree_optimization.jl")
    include(rhomb_folder * "observable_maps.jl")

    include(tbg_folder * "hf_presets.jl")
    include(tbg_folder * "hf_model.jl")
    include(tbg_folder * "hf_velocities.jl")
    include(tbg_folder * "rz_operator_ansatz.jl")
    include(tbg_folder * "optic_maps_correlated.jl")
    include(tbg_folder * "wrapper_lmc.jl")
    include(tbg_folder * "read_atlas.jl")
    
    include(common_folder * "separate_contributions.jl")
    include(proj_folder * "/cluster/cluster_tools.jl")
    include(proj_folder * "/cluster/merged_data.jl")

    export Params_rhombohedral, params_rhombohedral, xxx_lmc_presets, 
        lmc_presets, lmcshift_presets, lmcnoshift_presets, 
        xx_drude_presets, σxyahe_presets, qah_presets
    export Interpolated_params, Computation_params, Phase_diagram_params, Observable_computation
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
    export lmc_map, drude_map, ahe_map, map_eval, Ez_map_eval
    export script_path, slurm_submit_interpolations, slurm_submit_phasediagrams, slurm_submit_observable, data_merge, postprocessing
    export plot_phasediagrams, spinfull_plotbandsanddos, plot_cluster_bandsanddos
    export plot_drude, plot_ahe, plot_lmc, plot_lmcspin, aux_plot_obs
    export hf_valley_spin_hamiltonian, dhf_hamiltonian, rz, int_boundaries, paramsHF, ParamsHF
    export lmc_spin_sweep, lmc_orb_sweep, Self_consistent_data, drude_sweep
end