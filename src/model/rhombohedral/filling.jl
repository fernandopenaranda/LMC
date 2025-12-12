#= Remarks, the filling is tricky in the McCann-Koshino model for RH
the reason is that it is  k.p model and thus to determine the occupation
it requires in principle faithfull knowledge of the bandstructure at all k-points
in the BZ. However, provided the chemical potentials explored are small an in the 
order of magnitude of the Ez induced gap, then, we can integrate over the BZ safely
just by noting that for larger momenta than those giving $ϵ(k)> μ$ the bands in the
electrons (holes) grows (decreases) monotonously.
=#
"""  given a chemical potential value, the number of layers and the params presets
it computes the filling, it computes the filling over the k=bounds.
Important remark: if the k bounds are not the full BZ, the result will not correspond to the total
band filling. This is particularly important in continuum models.
In those cases we must ensure that the k bounds passed are larger than 
the k value k_μ where ϵ(k_μ) = μ. 
This correction will look as A_small/A_BZ * filling() 
where A_small is the area considered over the integration
"""
function rh_filling(N, p::Params_rhombohedral, μ; T = 0, ϵ = 1e-7, evals = 100)
    h(q) = abc_Nlayer(N, q, Params_rhombohedral(p, μ =0)) #unit_convention_two_packages_E
    cnst = p.γ1/(p.γ0 *√3/2)
    Δx = [-cnst-ϵ, cnst]   # this is a ratio of energies, convention independent
    Δy = [-cnst-ϵ, cnst] # the areas are unitless
    A_BZ = 8π^2/√3
    A_small = ((Δx[2]-Δx[1])*(Δy[2]-Δy[1]))
    return filling(h, μ, Δx, T, evals = evals) * A_small/A_BZ +
        (A_BZ-A_small)/A_BZ
end