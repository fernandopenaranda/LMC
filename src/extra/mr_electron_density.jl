"""
electron density in m^-2
"""

function electron_density(p, of::Matrix; evals = 100, T = 10)
    M, xmin, xmax = int_boundaries(p)
    integrand(q) = k_electron_density(p, q, of; T = T)
    val, err = Cubature.hcubature(integrand, [0, xmin[2]], [xmax[1]/2, xmax[2]];
        reltol=1e-5, abstol=0, maxevals = evals);
    return val / (4π^2*ang_to_m^2)
end

function k_electron_density(p, q, of; T = 10)
    h = hf_hamiltonian(p, q)
    ϵs, ψs = eigen(Matrix(h))
    return sum(f(ϵs, 0, T))
end


