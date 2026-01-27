function provide_folder_get_observables_atlas_reshuffled(folder, rs = "rs")
    str = "/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Projects/cluster_temp/exportdata"* "/" * string(folder)
       es_path = str * "/joined_es_"*rs*".csv"
       mus_path = str * "/joined_mus_"*rs*".csv"
       ns_path = str * "/joined_ns_"*rs*".csv"
       ofmats_path = str * "/joined_ofmats_"*rs*".csv"
       s = return_variables_reshuffled(es_path,mus_path, ns_path, ofmats_path )
       Self_consistent_data(s[1], s[2], s[3], s[4])
end


folder = 2281379#2281278
simportedvp = provide_folder_get_observables_atlas_reshuffled(folder,  "vp")
simportedqah = provide_folder_get_observables_atlas_reshuffled(folder,  "qah")


l = readparams(folder);
pimported = ParamsHF(paramsHF(1.05, 1,  1),
    μ = 0,
    nf = l.nf[1],
    λ = l.λ[1],
    M = l.M[1],
    v = l.v[1],
    vp = l.vp[1],
    γ =l.γ[1],
    sigmaz = l.sigmaz[1],
    sigmazlayerz = l.sigmazlayerz[1],
    U1 = l.U1[1],
    U2 = l.U2[1],
    J = 0*l.J[1],
    VP =l.VP[1],
    twovalleystwospins= l.twovalleystwospins[1]);

include("/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Projects/HeavyFermion_Optics/src/mr_module.jl")
include("/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Projects/HeavyFermion_Optics/src/mr_read_atlas")
include("/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Projects/HeavyFermion_Optics/src/mr_atlas_plotters.jl")

# # PID = "MR2604173-2281278"
# # PID = "MR2604409-2281278"
# # PID = "MR2605392-2281278"
# # PID = "MR2606874-2281278"
# # PID = "MR2607521-2281278"
# #PID = "MR2608158-2281278"
# # PID = "MR2604408-2281278"
# # PID = "MR2604411-2281278"
# # PID = "MR2606873-2281278"
# # PID = "MR2607520-2281278"
# # PID = "MR2608156-2281278" # T = 1 100000
# #  PID2 = "MR2608150-2281278"  # QAH 100000 1K  T = 1 100000
# PID = "MR2608157-2281278" # vp 100000 5K
# PID2= "MR2608158-2281278"# qah 100000 5K
# # PID ="MR2613437-2281278"
# # PID2 = "MR2613539-2281278"
# # No substrate
# PID = "MR2620964-2618196" #1K

folder = 2281379#2281278
simportedvp = provide_folder_get_observables_atlas_reshuffled(folder,  "vp")
simportedqah = provide_folder_get_observables_atlas_reshuffled(folder,  "qah")


ns, mcqah = mc_sweep(:x, :x, :x, pimported, simportedqah; 
    T = 1, τ = 200, evals = 5120, Ω_contr = true, omm_contr = true, fermi_surface = false)




    PID = string(ARGS[3])
rs = string(ARGS[4]) # vp or qah
evals = parse(Int, ARGS[5])
T = parse(Float64, ARGS[6])

job_id += 1
include("HF_Optics.jl")
include("degen_decision.jl") # to load the qah and vp structs
include("optic_maps_correlated.jl") # import the qah and vp structs
include("mr_module.jl")

p = read_struct_data(folder)
s = provide_folder_get_observables_atlas_reshuffled(folder, rs)


mc_map_correlated(folder,PID, job_id, :x,:x,:x, p, s.mus[job_id], diagm(s.ofmats[job_id]), s.ns[job_id],
    T = T, τ = 200, evals = evals, Ω_contr = true, omm_contr = true, rs = rs)


save_eval_presets(rs, evals, T, PID, job_id, folder) # save comp presets


function magneto_conductivity(i::Symbol, j::Symbol, k::Symbol, p, of; T = 2, τ = 100, evals = 100, Ω_contr = true, omm_contr = true, fermi_surface = false)
    τ *= femto_to_seconds
    M, xmin, xmax = int_boundaries(p)
    integrand(q) = k_linear_magnetorresistance(i, j, k, p, q, of; T = T, τ = τ, Ω_contr = Ω_contr, omm_contr = omm_contr, fermi_surface = fermi_surface)
    # integrand(q) = 1  # Volume of the BZ
    #integrand(q) = fermi_surface_integrand(p, q, of, p.μ, T) Fermi surface integral
    # val, err = Cubature.hcubature(integrand, [0, xmin[2]], [xmax[1]/2, xmax[2]]; # Right
    val, err = Cubature.hcubature(integrand, [0, xmin[2]], [xmax[1]/2, xmax[2]]; reltol=1e-5, abstol=0, maxevals = evals);
    return val
end

""" 
dhx and dhy do not include a momentum dependence
these have to change if we include a rz in the hamiltonian (for instance due to a sztz mass)
"""
function k_linear_magnetorresistance(i, j, k, p, q, of; T = 0, τ = 1e-15, Ω_contr = true, omm_contr = true, fermi_surface = false)
    err_dirs(i,j,k)
    h = hf_valley_spin_hamiltonian(p, q, of)
    ϵs, ψs = eigen(Matrix(h))
    dhx = dhf_hamiltonian(p, :x, of)
    dhy = dhf_hamiltonian(p, :y, of)
    rzmat = rz(ψs, p, q, of) * 3.3/2  # <-! # is the a0 needed !!! units and ξ???
    C = 1e3 * 2π * τ/(4π^2*ang_to_m^2)
    # the first 2π is the part of the reduced planck constant so σ_ijk is expressed in units of the quantum of conductance
    # 1e3 is needed to express the electron charge in meV so it cancels 1e3 Kb T in the denominator
    # 4π^2 comes from the k's diffs
    C * k_linear_mr_integrand(i, j, k, ϵs, ψs, rzmat, dhx, dhy, p.μ, T, Ω_contr = Ω_contr, omm_contr = omm_contr, fermi_surface = fermi_surface) 
end



#_________________________________________#_________________________________________
#_________________________________________#_________________________________________

function xxx_lmc_presets(p::ParamsHF, of; T = 10, τ = 200, evals = 100, ϵ = 1e-7,  
        berry_contribution = true, omm_contribution = true, fermi_surface = false, with_shift= false)
    unit_convention_two_packages_t = 1e-15
    h(q) = hf_valley_spin_hamiltonian(p, q, of)
    dhx(q) = dhf_hamiltonian(p, :x, of)    # q!!
    dhy(q) = dhf_hamiltonian(p, :y, of)
    dh(q) = [dhx(q), dhy(q)]
    dims = size(h([0.0,0.0]))
    dhxx(q) =  zeros(Float64, dims[1], dims[2]) # If there is no layer-odd sublattice potential v_xx is zero. the reason is that 
    # the only contribution to vxx may come from the non-local rz operator see our paper, but in our substrate model that
    # only breaks C2x it is zero, nevertheless even if other sublattice potentials are present its effect will be perturbative
    # as also demonstrated in our PRB.
    rzmat(q, ψs) = rz(ψs, p, q, of) * 3.3/2 
    
    M, xmin, xmax = int_boundaries(p)
    xbounds = [0, xmin[2]]
    ybounds = [xmax[1]/2, xmax[2]]
    # computation presets
    cpt = Transport_computation_presets(xbounds, ybounds, evals)
    # planar preset object (argument of `Optics_in_the_length_gauge.linear_magneto_conductivity_orbital`)
    planar_presets = Planar_σijk_presets_orbital(a0, :x,:x,:x, h, dh, dhxx, rz, τ*unit_convention_two_packages_t,
        T, cpt, berry_contribution, omm_contribution, fermi_surface, with_shift)
    return planar_presets
end
