# BANDS
function abck_Nbands(N, num_points, p)
    k_vecs = abcbz_path_gamma_k_m_gamma(num_points, p)
    # Rs = ([√3/2,-1/2], [0,1.0])
    # sgnum = 17                
    # k_vecs = kpath(Rs, sgnum, num_points, [:Γ, :M, :K, :Γ])
    es = zeros(ComplexF64, 2,length(k_vecs))
    for i in 1:length(k_vecs)
        e = eigen(Matrix(abc_Nlayer(N, [k_vecs[i][1], k_vecs[i][2]], p)))
        es[:,i] = e.values 
    end
    es
end

function abck_bands(num_points, p)
    k_vecs = abcbz_path_gamma_k_m_gamma(num_points, p)
    es = zeros(ComplexF64, 2,length(k_vecs))
    for i in 1:length(k_vecs)
        e = eigen(Matrix(abc_pentalayer([k_vecs[i][1], k_vecs[i][2]], p)))
        es[:,i] = e.values 
    end
    es
end
#=_____________________________________________________________________________________________
Kpaths
____________________________________________________________________________________________=#
"""
path very close to the K_ξ point, where the two band model works.
momenta is adimensional
"""
function abcbz_path_gamma_k_m_gamma(nk::Int, p)
    coef = p.γ1/p.γ0 *√3/2
    function interpolate(p1, p2, n)
        [(p1[1] + (p2[1] - p1[1]) * t, p1[2] + (p2[2] - p1[2]) * t) for t in range(0, stop=1, length=n+1)[1:end-1]]
    end
    path = interpolate((-coef, 0), (coef, 0), 2nk)
    return path
end