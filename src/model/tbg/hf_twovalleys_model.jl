"""
Kramers intervalley coherent groundstate at nu =0
We take the one-shot approach. 
cc -> J σyτy
ff -> U1 σyτy
"""
# function hf_kivc_ν0(p, k)
#     h = hf_2vhamiltonian(ParamsHF(p, nf = 1, n = 0, VP = false), k)
#     kivc!(h, p)
#     return h
# end

function kivc!(h_mat, p)
    dim = size(h_mat, 1) ÷ 2
    # density density local (no hartree) Note that it is finite because at charge neutrality one valley is above and the other is below.
    h_mat[1:2, dim+1:dim+2] .+= -1im * p.U1/2 .* σy #\sigmay \tauy
    h_mat[dim+1:dim+2, 1:2] .+= 1im * p.U1/2 .* σy
    # exchange interaction in the conduction bands corresponding to the gamma_1+gamma_2 sector
    g = zeros(Float64, 2)
    ham_dim = size(h_mat, 1)
    g1, g2 = bravais_vectors(p)
    cc_mat = 1im * p.J .* σy
    for i in 3:4:dim-3
        h_mat[i+2:i+3, dim+i+2:dim+i+3] .+= cc_mat
        h_mat[dim+i+2:dim+i+3,i+2:i+3] .+= -cc_mat
    end
end


# hf_plotbands(ParamsHF(p, c2ybreakingmass = 0, sigmaz = 0,  μ = 0, vafek = 0, J = 0, U1 = 35, VP = false, twovalleys = true, KIVC = true))