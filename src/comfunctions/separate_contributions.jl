"""
function that computes the k_resolved OMM and/or Berry curvature contributions to the LMC
"""
function klmc(pR::Params_rhombohedral, i,j,k, h, dh, ddh, rz, τ, T; 
        Ω_contr = true, omm_contr = true, fermi_surface = false, 
        with_shift = true, points = 10)
    integrand(q) = Optics_in_the_length_gauge.k_linear_magneto_conductivity_orbital(i, j, k, h, dh, ddh, rz, q; T = T, τ = τ, 
        Ω_contr = Ω_contr, omm_contr = omm_contr, fermi_surface = fermi_surface, with_shift = with_shift)
    return evalmat(pR, integrand, points)
end

"""
k -resolved Ωz
"""
function kresolved_Ωz(pR::Params_rhombohedral, p::Planar_σijk_presets_orbital; points = 10)
    integrand(q) = k_Omegaz(p, q)
    return evalmat(pR, integrand, points)
end

"""
function k-resolved Ωin
"""
function kresolved_Ωin(pR::Params_rhombohedral, p::Planar_σijk_presets_orbital; points = 10)
    integrand(q) = k_Omegain(p, q)
    return evalmat(pR, integrand, points)
end

function kresolved_dOMM(pR::Params_rhombohedral, p::Planar_σijk_presets_orbital; points = 10) 
    integrand(q) = k_d_OMM(p, q)
    return evalmat(pR, integrand, points)
end

"""evaluator of a function f in the BZ"""
function evalmat(p, f, points)
    cnst = p.γ1/(p.γ0 *√3/2)
    kxs = range(-cnst ,cnst, length = points)
    kys = range(-cnst, cnst, length = points)
    Zs = [f([kx, ky]) for kx in kxs, ky in kys]
    return kxs, kys, Zs
end