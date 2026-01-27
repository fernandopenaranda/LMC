"""
    wrapper to Optics_in_the_length_gauge functions for a LMC conductivity calculations
    in the THFM.
    If there is no layer-odd sublattice potential v_xx is zero. the reason is that 
    the only contribution to vxx may come from the non-local repr of the rz operator 
    see our paper, but in our substrate model that only breaks C2x it is zero
    nevertheless even if other sublattice potentials are present its effect will be 
    perturbative as also demonstrated in our PRB.
"""
function xxx_lmc_presets(p::ParamsHF, of; T = 1, τ = 200, evals = 100, ϵ = 1e-7,  
        berry_contribution = true, omm_contribution = true, fermi_surface = false,
        with_shift= false)
    unit_convention_two_packages_t = 1e-15
    h(q) = hf_valley_spin_hamiltonian(p, q, of)
    dhx(q) = dhf_hamiltonian(p, :x, of)    # q!!
    dhy(q) = dhf_hamiltonian(p, :y, of)
    dh(q) = [dhx(q), dhy(q)]
    dims = size(h([0.0,0.0]))
    dhxx(q) =  zeros(Float64, dims[1], dims[2]) # 
    rzmat(q, ψs) = rz(ψs, p, q, of) * 3.3/2 

    M, xmin, xmax = LMC.int_boundaries(p)
    xbounds = [0, xmin[2]]
    ybounds = [xmax[1]/2, xmax[2]]
    # computation presets
    cpt = Transport_computation_presets(xbounds, ybounds, evals)
    # planar preset object (argument of `Optics_in_the_length_gauge.linear_magneto_conductivity_orbital`)
    planar_presets = Planar_σijk_presets_orbital(a0, :x,:x,:x, h, dh, dhxx, rzmat, τ*unit_convention_two_packages_t,
        T, cpt, berry_contribution, omm_contribution, fermi_surface, with_shift)
    return planar_presets
end



function xxx_lmc_spin_presets(p::ParamsHF, of; T = 1, τ = 200, evals = 100, ϵ = 1e-7,  
    berry_contribution = true, omm_contribution = true, fermi_surface = false,
    with_shift= false)
unit_convention_two_packages_t = 1e-15
h(q) = hf_valley_spin_hamiltonian(p, q, of)
dhx(q) = dhf_hamiltonian(p, :x, of)    # q!!
dhy(q) = dhf_hamiltonian(p, :y, of)
dh(q) = [dhx(q), dhy(q)]
dims = size(h([0.0,0.0]))
dhxx(q) =  zeros(Float64, dims[1], dims[2]) # 
rzmat(q, ψs) = rz(ψs, p, q, of) * 3.3/2 

M, xmin, xmax = LMC.int_boundaries(p)
xbounds = [0, xmin[2]]
ybounds = [xmax[1]/2, xmax[2]]
# computation presets
cpt = Transport_computation_presets(xbounds, ybounds, evals)
# planar preset object (argument of `Optics_in_the_length_gauge.linear_magneto_conductivity_orbital`)
planar_presets = Planar_σijk_presets_orbital(a0, :x,:x,:x, h, dh, dhxx, rzmat, τ*unit_convention_two_packages_t,
    T, cpt, berry_contribution, omm_contribution, fermi_surface, with_shift)
return planar_presets
end

"""
semiclassical expression for the field independent drude conductivity
Units: [e^2/h]
returns the presets struct that set the 
    `Optics_in_the_length_gauge.drude_conductivity` calculation
"""
function xx_drude_presets(p::ParamsHF, of; T = 1, τ = 200, evals = 100, ϵ = 1e-7)
        unit_convention_two_packages_t = 1e-15
    h(q) = hf_valley_spin_hamiltonian(p, q, of)
    dhx(q) = dhf_hamiltonian(p, :x, of)  
    M, xmin, xmax = LMC.int_boundaries(p)
    xbounds = [0, xmin[2]]
    ybounds = [xmax[1]/2, xmax[2]]
    cpt = Transport_computation_presets(xbounds, ybounds, evals)
    return Drude_presets(a0, :x, :x, h, dhx, T, 
        τ*unit_convention_two_packages_t, cpt)
end
