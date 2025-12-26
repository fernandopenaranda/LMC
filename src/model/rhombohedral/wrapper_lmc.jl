"""
semiclassical expression for the linear in B magnetoconductivity in the presence of constant 
electric and magnetic field and electric fields in the plane.

The term vij is zero if there is no substrate. 
If there is a dependency on rz which depends on sin(k) will give rise to a vij ≠ 0

Units: [e^2/h * 1/T]
τ is the scattering time in seconds
returns the presets struct that set the `Optics_in_the_length_gauge.linear_magneto_conductivity_orbital` calculation
The model is written in meV, fs, Angstroms, and K.
The module Optics_in_the_length_gauge is written in eV, s, K, so there is a unit 
    convention fixer
"""
xxx_lmc_presets(N, μ, p::Params_rhombohedral; kws...) = 
    xxx_lmc_presets(N, Params_rhombohedral(p, μ = μ); kws...)
xxx_lmc_presets(N, μ, ξ, p::Params_rhombohedral; kws...) = 
    xxx_lmc_presets(N, Params_rhombohedral(p, ξ = ξ, μ = μ); kws...)
function xxx_lmc_presets(N, p::Params_rhombohedral; T = 10, τ = 200, evals = 100, ϵ = 1e-7,  
        berry_contribution = true, omm_contribution = true, fermi_surface = false, with_shift= false)
    unit_convention_two_packages_E = 1e-3
    unit_convention_two_packages_t = 1e-15
    h(q) = abc_Nlayer(N, q, p) .* unit_convention_two_packages_E
    dhx(q) = dhxNlg(N, q, p)   .* unit_convention_two_packages_E
    dhy(q) = dhyNlg(N, q, p)   .* unit_convention_two_packages_E
    dh(q) = [dhx(q), dhy(q)]
    dhxx(q) = dhxxNlg(N, q, p) .* unit_convention_two_packages_E
    
    rz(q, ψs) = rzNlg(N, ψs)  
    # integral bounds (around the valley)
    cnst = p.γ1/(p.γ0 *√3/2)
    xbounds = [-cnst-ϵ, cnst]   # this is a ratio of energies, convention independent
    ybounds = [-cnst-ϵ, cnst]   
    # computation presets
    cpt = Transport_computation_presets(xbounds, ybounds, evals)
    # planar preset object (argument of `Optics_in_the_length_gauge.linear_magneto_conductivity_orbital`)
    planar_presets = Planar_σijk_presets_orbital(a0, :x,:x,:x, h, dh, dhxx, rz, τ*unit_convention_two_packages_t,
        T, cpt, berry_contribution, omm_contribution, fermi_surface, with_shift)
    return planar_presets
end

"""
semiclassical expression for the field independent drude conductivity

Units: [e^2/h]
returns the presets struct that set the `Optics_in_the_length_gauge.drude_conductivity` calculation
"""
# drude_conductivity(p::Drude_presets) = drude_conductivity(p.h, p.dhi, p.T, p.τ, 
#         p.computation.xbounds, p.computation.ybounds, p.computation.evals)

xx_drude_presets(N, μ, p::Params_rhombohedral; kws...) = 
    xx_drude_presets(N, Params_rhombohedral(p, μ = μ); kws...)
xx_drude_presets(N, μ, ξ, p::Params_rhombohedral; kws...) = 
    xx_drude_presets(N, Params_rhombohedral(p, ξ = ξ, μ = μ); kws...)
function xx_drude_presets(N, p::Params_rhombohedral; 
        T = 10, τ = 200, evals = 100, ϵ = 1e-7)
    unit_convention_two_packages_E = 1e-3
    unit_convention_two_packages_t = 1e-15
    h(q) = abc_Nlayer(N, q, p) .* unit_convention_two_packages_E
    dhx(q) = dhxNlg(N, q, p)   .* unit_convention_two_packages_E
    # integral bounds (around the valley)
    cnst = p.γ1/(p.γ0 * √3/2)
    xbounds = [-cnst-ϵ, cnst]   # this is a ratio of energies, 
    # convention independent
    ybounds = [-cnst-ϵ, cnst]   
    # computation presets
    cpt = Transport_computation_presets(xbounds, ybounds, evals)
    # planar preset object (argument of 
    # `Optics_in_the_length_gauge.linear_magneto_conductivity_orbital`)
    return Drude_presets(a0, :x,:x, h, dhx, T, 
        τ*unit_convention_two_packages_t, cpt)
end

#_______________________________________________________________________________________
# More wrappers
#_______________________________________________________________________________________
lmc_presets(N, μ,ξ, p::Params_rhombohedral; kws...) = 
     xxx_lmc_presets(N, μ, ξ, p; kws...)
lmcnoshift_presets(μ,ξ, evals) = 
    xxx_lmc_presets(N, μ, ξ, p; evals = evals, T = T, 
    berry_contribution = true, omm_contribution = true, fermi_surface = false, 
    with_shift = false)
lmcshift_presets(μ,ξ) = 
    xxx_lmc_presets(N, μ, ξ, p; evals = evals, T = T, fermi_surface = false, 
    with_shift = true)

#_______________________________________________________________________________________
# DOS
#_______________________________________________________________________________________
""" density of states call to Optics_in_the_length_gauge. Arguments in meV """
c_dos(p::Planar_σijk_presets_orbital, μ::Number; η = 0.005, evals = 100) = 
    c_dos(p::Planar_σijk_presets_orbital, [μ]; η = η, evals = evals)

c_dos(p::Planar_σijk_presets_orbital, μlist::Array; η = 0.005, evals = 10000) = 
    Optics_in_the_length_gauge.dos(p.a0, p.h, p.computation.xbounds, p.computation.ybounds, μlist ./ 1e3, η = η/1e3, evals = evals)
#_______________________________________________________________________________________
# QAH
#_______________________________________________________________________________________
""" quantum anomalous Hall response presets builder """
σxyahe_presets(N, μ, ξ, p::Params_rhombohedral; kws...) = qah_presets(N, :x, :y, μ, ξ, p; kws...)
qah_presets(N, dirJ, dirE, μ, ξ, p::Params_rhombohedral; kws...) = qah_presets(N, dirJ, dirE, Params_rhombohedral(p, ξ = ξ, μ = μ); kws...)
function qah_presets(N, dirJ, dirE, p::Params_rhombohedral; T = 10, evals = 100, ϵ = 1e-7) 
    unit_convention_two_packages_E = 1e-3
    unit_convention_two_packages_t = 1e-15
    h(q) = abc_Nlayer(N, q, p) .* unit_convention_two_packages_E
    dhx(q) = dhxNlg(N, q, p)   .* unit_convention_two_packages_E
    dhy(q) = dhyNlg(N, q, p)   .* unit_convention_two_packages_E
    dh(q) = [dhx(q), dhy(q)]
    rz(q, ψs) = rzNlg(N, ψs)  
    # integral bounds (around the valley)
    cnst = p.γ1/(p.γ0 *√3/2)
    xbounds = [-cnst-ϵ, cnst]   # this is a ratio of energies, convention independent
    ybounds = [-cnst-ϵ, cnst]   
    # computation presets
    cpt = Transport_computation_presets(xbounds, ybounds, evals)
    return AH_presets(a0, dirJ, dirE, h, dh, T, cpt)
end
#_________________________________________________________________________________________
# k - resolved calculations
#_________________________________________________________________________________________
"""kresolved calcuations"""
kresolvedlmc(pR::Params_rhombohedral, p::Planar_σijk_presets_orbital; Ω_contr = true, omm_contr = true, fermi_surface = false, with_shift = true, points = 10) =
    klmc(pR, p.dirJ,p.dirE,p.dirB, p.h, p.nabla_h, p.nabla_nabla_h, p.rz, p.τ, p.T, 
        Ω_contr = Ω_contr, omm_contr = omm_contr, fermi_surface = fermi_surface, with_shift = with_shift, points = points)
"""kresolved inplane berry curvature"""
k_Omegain(p::Planar_σijk_presets_orbital, q) = Optics_in_the_length_gauge.k_Ωi_fs(p.dirJ, p.dirE, p.h, p.nabla_h, p.rz, q, p.T)
"""kresolved outofplane berry curvature"""
k_Omegaz(p::Planar_σijk_presets_orbital, q) = Optics_in_the_length_gauge.k_Ωxy_fn(p.dirJ, p.dirE, p.h, p.nabla_h, p.rz, q, p.T)
"""kresolved omm"""
k_d_OMM(p::Planar_σijk_presets_orbital, q) = Optics_in_the_length_gauge.k_d_OMM_fs(p.dirJ, p.dirE, p.h, p.nabla_h, p.rz, q, p.T)


#____________________________________________________________________________________________

# function magneto_conductivity(N, i::Symbol, j::Symbol, k::Symbol, p::Params_rhombohedral; 
#         T = 2, τ = 200, evals = 100, Ω_contr = true, omm_contr = true, fermi_surface = false)
#     τ *= femto_to_seconds
#     cnst = 1*p.γ1/p.γ0 *√3/2 #*20
#     cnst *= 1
#     xmin, xmax = [-cnst, -cnst],[cnst, cnst]
#     integrand(q) = Nk_linear_magnetorresistance(N, i, j, k, p, q; T = T, τ = τ, 
#         Ω_contr = Ω_contr, omm_contr = omm_contr, fermi_surface = fermi_surface)
#     val, err = Cubature.hcubature(integrand, [-cnst+cnst/1e3,-cnst], [cnst, cnst];
#          reltol = 1e-5, abstol=0, maxevals = evals);
#     return val
# end
# [xbounds[1], ybounds[1]], [xbounds[2], ybounds[2]]
# [-cnst+cnst/1e3,-cnst], [cnst, cnst]
# """
# semiclassical expression for the linear in B magnetoconductivity in the presence of constant 
# electric and magnetic field and electric fields in the plane.
# Units: [e^2/h] Already multiplied by B. 
# B in Teslas
# """
# function magneto_conductivity_withshift(N, i::Symbol, j::Symbol, k::Symbol, p::Params_rhombohedral;
#          B =10, T = 2, τ = 200, evals = 100, Ω_contr = true, omm_contr = true, fermi_surface = false)
#     τ *= femto_to_seconds
#     cnst = 1*p.γ1/p.γ0 *√3/2
#     cnst *= 1
#     xmin, xmax = [-cnst, -cnst],[cnst, cnst]
#     integrand(q) = Nk_linear_magnetorresistance_withshift(N, i, j, k, p, q, B; T = T, 
#         τ = τ, Ω_contr = Ω_contr, omm_contr = omm_contr, fermi_surface = fermi_surface)
   
#     fdim = 4
#     val, err = Cubature.hcubature(fdim, (q,v) -> v[:] = integrand(q), [-cnst+cnst/1e3,-cnst], [cnst, cnst]; 
#         reltol = 1e-5, abstol=0, maxevals = evals);
#     return B*val[1] + 0*(val[4]/val[3])*val[2] #[σ, vij, N0, δμ]
# end

#____________________________________________________________________________________________



# """
#     k_linear_magnetorresistance
# integrand of σ_ijk with all the prefactors
# About units:
#  the first 2π is the part of the reduced planck constant so σ_ijk is expressed in units of the quantum of conductance
#  1e3 is needed to express the electron charge in meV so it cancels 1e3 Kb T in the denominator
#  4π^2 comes from the k's diffs
# """
# function k_linear_magnetorresistance(i::Symbol, j::Symbol, k::Symbol, p::Params_rhombohedral, q; 
#         T = 2, τ = 1e-15, Ω_contr = true, omm_contr = true, fermi_surface = false)
#     h = abc_pentalayer(q, p)
#     ϵs, ψs = eigen(Matrix(h))
#     dhx = dhx5lg(q, p)
#     dhy = dhy5lg(q, p)
#     dhxx = dhxxNlg(5, q, p)
#     rzmat = rz5lg(ψs) # el 4 viene de N-1 layers en la distancia entre layers 
#     C = 1e3 * 2π * τ/ (4π^2*ang_to_m^2)
#     C * k_linear_mr_integrand(p, i, j, k, ϵs, ψs, rzmat, dhx, dhy, dhxx, 0, T, 
#         Ω_contr = Ω_contr, omm_contr = omm_contr, fermi_surface = fermi_surface)
# end

# """
#     k_linear_magnetorresistance
# integrand of σ_ijk with all the prefactors
# About units:
#  the first 2π is the part of the reduced planck constant so σ_ijk is expressed in units of the quantum of conductance
#  1e3 is needed to express the electron charge in meV so it cancels 1e3 Kb T in the denominator
#  4π^2 comes from the k's diffs
# """
# function Nk_linear_magnetorresistance(N, i::Symbol, j::Symbol, k::Symbol, p::Params_rhombohedral, q; 
#         T = 2, τ = 1e-15, Ω_contr = true, omm_contr = true, fermi_surface = false)
#     h = abc_Nlayer(N, q, p)
#     ϵs, ψs = eigen(Matrix(h))
#     dhx = dhxNlg(N, q, p)
#     dhy = dhyNlg(N, q, p)
#     dhxx = dhxxNlg(N, q, p)
#     rzmat =rzNlg(N, ψs) # el 4 viene de N-1 layers en la distancia entre layers 
#     C = 1e3 * 2π * τ/ (4π^2*ang_to_m^2)
#     C * k_linear_mr_integrand(p, i, j, k, ϵs, ψs, rzmat, dhx, dhy, dhxx, 0, T, 
#         Ω_contr = Ω_contr, omm_contr = omm_contr, fermi_surface = fermi_surface)
# end

# function Nk_linear_magnetorresistance_withshift(N, i::Symbol, j::Symbol, k::Symbol, p::Params_rhombohedral, q, B; 
#         T = 2, τ = 1e-15, Ω_contr = true, omm_contr = true, fermi_surface = false)
#     h = abc_Nlayer(N, q, p)
#     ϵs, ψs = eigen(Matrix(h))
#     dhx = dhxNlg(N, q, p)
#     dhy = dhyNlg(N, q, p)
#     dhxx = dhxxNlg(N, q, p)
#     vxx = vel(ψs, dhxx) * ang_to_m^2/ ħ_mev_s
#     ry = r(ϵs, ψs, dhy) * ang_to_m
#     vy = vel(ψs, dhy) * ang_to_m/ ħ_mev_s
#     rzmat =rzNlg(N, ψs) # el 4 viene de N-1 layers en la distancia entre layers 
#     C = 1e3 * 2π * τ
#     dif_k_units = 1/(4π^2*ang_to_m^2) 
    
#     SA[dif_k_units * C * k_linear_mr_integrand(p, i, j, k, ϵs, ψs, rzmat, dhx, dhy, dhxx, 0, T, 
#         Ω_contr = Ω_contr, omm_contr = omm_contr, fermi_surface = fermi_surface),
#         C* dif_k_units .* vij_shift(ϵs, T, vxx),
#         k_dos(p, ϵs, T), #units of E * m^2
#         δμ_shift(i, B, ϵs, T, vy, ry, rzmat * ang_to_m)] # 1/(m^2 e)
# end

# "the term with vij is only valid for sigma xxx"
# function k_linear_mr_integrand(p::Params_rhombohedral, i, j, k, ϵs, ψs, rzmat, dhx, dhy, dhxx, μ, T;
#          Ω_contr = true, omm_contr = true, fermi_surface = false) # units meters, meV, seconds
#     omega = Ω(ϵs)
#     Δx = Δ(ψs, dhx) * ang_to_m
#     Δy = Δ(ψs, dhy) * ang_to_m
#     rx = r(ϵs, ψs, dhx) * ang_to_m
#     ry = r(ϵs, ψs, dhy) * ang_to_m
#     vx = vel(ψs, dhx) * ang_to_m/ ħ_mev_s
#     vy = vel(ψs, dhy) * ang_to_m/ ħ_mev_s
#     vxx = vel(ψs, dhxx) * ang_to_m^2/ ħ_mev_s
#     rzmat *= ang_to_m
#     Ω_switch = ifelse(Ω_contr == true, 1, 0)
#     omm_switch = ifelse(omm_contr == true, 1, 0)
#     if fermi_surface == true
#         return sum(d_f(ϵs, 0, T))
#     else
#         return sum(d_f(ϵs, 0, T) .*
#             (omm_switch .* mr_omm(i, j, omega, rx, ry, vx, vy, Δx, Δy, rzmat) + 
#             Ω_switch .* mr_Ω(i, j, k, rzmat, rx, ry, vx, vy) + 
#             - mr_vij(i, vy, rzmat, vxx) #only valid in the xx direction
#             ))
#     end
# end


# """"
# correction due to switching into the canonical ensemble
# """
# vij_shift(ϵs, T, vij) = sum(d_f(ϵs, 0, T) .* real(diag(vij)))
# mr_vij(i, vj,rz, vij) = real(OMM(i, vj, rz) .* diag(vij))

# δμ_shift(i, B, ϵs, T, vj, rj, rz) = sum(d_f(ϵs, 0, T) .* (B .* OMM(i, vj, rz)) + 
#    (B .* Ωin(i, rj, rz)) .* f(ϵs, 0, T)/ħ_mev_s) #units 1/e m^2

# """
#     σxx
# """
# function in_plane_bindependent_conductivity(p::Params_rhombohedral; T = 2, τ = 200, evals = 100)
#     τ *= femto_to_seconds
#     cnst = 1*p.γ1/p.γ0 *√3/2 #*20
#     cnst *= 1
#     xmin, xmax = [-cnst, -cnst],[cnst, cnst]
#     integrand(q) = k_in_plane_bindependent_conductivity(p, q; T = T, τ = τ)
#     val, err = Cubature.hcubature(integrand, [-cnst+cnst/1e3,-cnst], [cnst, cnst]; 
#         reltol=1e-5, abstol=0, maxevals = evals);
#     return val
# end

# function in_plane_bindependent_conductivity(N, p::Params_rhombohedral; T = 2, τ = 200, evals = 100)
#     τ *= femto_to_seconds
#     cnst = 1*p.γ1/p.γ0 *√3/2 #*20
#     cnst *= 1
#     xmin, xmax = [-cnst, -cnst],[cnst, cnst]
#     integrand(q) = k_in_plane_bindependent_conductivity(N, p, q; T = T, τ = τ)
#     val, err = Cubature.hcubature(integrand, [-cnst+cnst/1e3,-cnst], [cnst, cnst]; 
#         reltol=1e-5, abstol=0, maxevals = evals);
#     return val
# end

# """
#     k_linear_magnetorresistance
# integrand of σ_xx with all the prefactors
# """
# function k_in_plane_bindependent_conductivity(p::Params_rhombohedral, q; T = 2, τ = 1e-15)    
#     h = abc_pentalayer(q, p)
#     ϵs, ψs = eigen(Matrix(h))
#     dhx = dhx5lg(q, p)
#     C = 2π * τ* ħ_mev_s/(4π^2*ang_to_m^2) 
#     # 1e3 is needed to express the electron charge in meV so it cancels 1e3 Kb T in the denominator
#     # 4π^2 comes from the k's diffs
#     C * k_in_plane_bindependent_conductivity_integrand(:x, ϵs, ψs, dhx, p.μ, T)
# end

# """
#     k_linear_magnetorresistance
# integrand of σ_xx with all the prefactors
# """
# function k_in_plane_bindependent_conductivity(N, p::Params_rhombohedral, q; T = 2, τ = 1e-15)  
#     h = abc_Nlayer(N,q, p)
#     ϵs, ψs = eigen(Matrix(h))
#     dhx = dhxNlg(N,q, p)
#     C = 2π * τ* ħ_mev_s/(4π^2*ang_to_m^2) 
#     # 1e3 is needed to express the electron charge in meV so it cancels 1e3 Kb T in the denominator
#     # 4π^2 comes from the k's diffs
#     C * k_in_plane_bindependent_conductivity_integrand(:x, ϵs, ψs, dhx, p.μ, T)
# end

