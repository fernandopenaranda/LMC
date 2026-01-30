function bands(p::ParamsHF, num_points; kws...)
    k_vecs = k_path(p, num_points)
    # fdim = p.twovalleys == false ? 1 : 2
    es = zeros(ComplexF64, Int(ham_matrix_size(p)),length(k_vecs))
    whichham = p.twovalleys == false ? hf_hamiltonian : 
    for i in 1:length(k_vecs)
        e = eigen(Matrix(whichham(p, [k_vecs[i][1], k_vecs[i][2]])))
        es[:,i] = e.values
    end
    es
end

bands(p::ParamsHF, of, num_points::Int64; kws...) = bands(p, k_path(p, num_points), of; kws...)

bands(p::ParamsHF, k_vecs, of::Vector; kws...) = bands(p, k_vecs, diagm(of); kws...)

function bands(p::ParamsHF, k_vecs, of::Matrix; kws...)
    # fdim = p.twovalleys == false ? 1 : 2
    es = zeros(ComplexF64, Int(ham_matrix_size(p)), length(k_vecs))
    if p.twovalleys == true
        whichham = hf_twovalleyshamiltonian
    elseif p.twovalleystwospins == true
        whichham = hf_valley_spin_hamiltonian
    end
    for i in 1:length(k_vecs)
        e = eigen(Matrix( whichham(p, [k_vecs[i][1], k_vecs[i][2]], of)))
        es[:,i] = e.values
    end
    es
end


function hartree_eigs(p, k; kws...)
    whichham = p.twovalleys == false ? hf_hamiltonian : hf_twovalleyshamiltonian
    l = eigen(Matrix(whichham(p, [k[1], k[2]]; kws...)))
    inds = sortperm(l.values, by=real)
    return l.values[inds], l.vectors[:, inds]
end


function sc_eigs(p, k, of; kws...)
    return eigen!(Matrix(hf_valley_spin_hamiltonian(p, [k[1], k[2]], of; kws...)))
end

function sc_eigvals!(evals, p, k, of; kws...)
    evals .= sort(real.(LinearAlgebra.LAPACK.geev!('N', 'N', 
        Matrix(hf_valley_spin_hamiltonian(p, [k[1], k[2]], of; kws...)))[1]))
end 

function sc_eigs!(evals, evecs, p, k, of; kws...)
    F = eigen!(Matrix(hf_valley_spin_hamiltonian(p, [k[1], k[2]], of; kws...)))
    evals .= F.values
    evecs .= F.vectors
end

"""
sparse diagonalization
"""
function sc_eigs_flatbands(p, k, of; kws...)
    whicham = (p.twovalleys == true ? hf_twovalleyshamiltonian : hf_valley_spin_hamiltonian)
    l = eigs(whicham(p, [k[1], k[2]], of; kws...), )
end

"""
returns the number of f electrons below filling n.
Note that n is made out of c and f electrons. Since the f electrons
are shared among dispersive and flat bands at n = 4, n_f<2
"""
function fdensity(p; kpoints = 42, kws...)
    kpoints = sixthMBZ(p; kpoints = kpoints )  # 1/6th of the MBZ
    hdim = ham_matrix_size(p)
    data = zeros(Float64, 2length(kpoints)); 
    for (i, k) in enumerate(kpoints)
        eigs = real.(hartree_eigs(p, k; kws...)[1])
        data[2i-1:2i] .= eigs[hdim÷2:hdim÷2+1] # flat band states
    end
    ind = Int(round(p.n/4 * length(kpoints)))
    auxmu = sort(data)
    mu = auxmu[ifelse(ind == -length(kpoints), 1, length(kpoints) + ind)]
    println("μ: ", mu," [meV]")
    nf = 0.0
    meshdim = length(kpoints)
    for k in kpoints
        nf += fdensityink(p, k, mu, meshdim; kws...)
    end
    return nf
end

function fdensityink(p, k, mu, meshdim; kws...) # works for a single and a two valley model
    whichham = p.twovalleys == false ? hf_hamiltonian : hf_twovalleyshamiltonian
    l = eigen(Matrix(whichham(p, [k[1], k[2]]; kws...)))
    meshsizerenorm = (6 * meshdim - 5)/ifelse(k == [0.,0.], 1, 6)  
    inds = sortperm(l.values, by=real)
    mat = l.vectors[:, inds]
    
    aux = spzeros(ComplexF64, size(mat,1), size(mat,1))
    fprojector = spzeros(ComplexF64, size(mat,1), size(mat,1))
    fprojector[1:2,1:2]  .= σ0
    halfdim = size(mat,1)÷2
    if p.twovalleys == true 
        fprojector[1+halfdim:2+halfdim,1+halfdim:2+halfdim] .+= σ0
    else nothing end
    for i in 1:size(mat,1)
        if l.values[inds][i] ≤ mu
            aux .+= (fprojector * mat[:,i]) * (fprojector * mat[:,i])'
        else nothing end
    end
    return sum(diag(aux))/ meshsizerenorm 
    # return diag(aux)/ meshsizerenorm 
end
"""
    computes O^f matrix. See HF appendix.
    when choosing a GS we will force this matrix to be of a particular type:
        if VP then block diagonal in the valleys
    Of(ParamsHF(p, sigmaz=0,vafek= 0, J = -8, U1 = 35, VP = true, twovalleys = true), -0)
    # return hf_plotbands(ParamsHF(p, μ = mu))
"""
Of(p::ParamsHF, mu, of; kws...) = expected_value_op(p, mu, of, 1I; kws...)

function expected_value_op(p::ParamsHF, mu, of, op; evals = 800)
    g1, g2 = bravais_vectors(p)
    κ1 = κ(g1, g2)
    ofmat = zeros(Float64, 8)
    halfdim = ham_matrix_size(p) ÷ 2
    projs = SA[1,2,1+halfdim÷2,2+halfdim÷2,1+halfdim,2+halfdim,1+3halfdim÷2,2+3halfdim÷2]
    vals = zeros(Float64, 2halfdim)
    vecs = zeros(ComplexF64, 2halfdim, 2halfdim)
    integrand(k) = k_expected_value_op!(vals, vecs, ofmat, projs, p, k, mu, of, op)
    
    M, xmin, xmax = int_boundaries(p)
    Δx = [0, xmin[2]]
    Δy = [xmax[1]/2, xmax[2]]
    
    val, err = Cubature.hcubature(8, (x,v) -> v[:] = integrand(x), Δx, Δy; maxevals = evals)
    return val /( (Δy[2] - Δx[2]) * (Δy[1] - Δx[1]))
end

function k_expected_value_op!(vals, vecs, ofmat, projs, p, k, mu, of, op; kws...)
    sc_eigs!(vals, vecs, p, [k[1], k[2]], of; kws...)
    # vals, vecs = eigen!(Matrix(hf_valley_spin_hamiltonian(p, [k[1], k[2]], of; kws...)))
    mat = vecs[:, findall(es-> es < mu, vals)]
    s = zero(eltype(ofmat)) 
    for (i1, p1) in enumerate(projs)
        for (i2, p2) in enumerate(projs)
            s *= 0.0
            op_val = op[i1,i2] 
            s += mat[p1,:]' * op_val * (mat[p2,:])
            if i1 == i2
                ofmat[i1] = real(s)
            else nothing end
        end
    end
    return ofmat
end

# Bandstructure
function k_path(p::ParamsHF, num_points; ϵ = 0)
    g1, g2 = bravais_vectors(p)
    K1 =  κ(g1, g2)
    K2 =  κ(g2, g1)
    # Initialize an array to store the line points
    # line_points = Array{Tuple{Float64, Float64}}(undef, num_points+1)
    # Compute the line points and store them in the array
    line_points = []
    for i in 0:(√3/2+1+1/2) * num_points+1
        if i ≤ num_points
        x = K1[1] - i * (1+ϵ) * K1[1]/num_points
        y = K1[2] + i *(1+ϵ) * K2[2]/num_points

        elseif i > num_points && i ≤ (√3/2+1) * num_points
        x = 0.0 + (i-num_points) * (1+ϵ) * K1[1]/(√3/2*num_points)
        y = 0.0
        else
        x = K1[1]
        y = 0 + (i-(√3/2+1) * num_points) *(1+ϵ) * K2[2]/(num_points/2)
        # line_points[i+1] = (x, y)
        end
    push!(line_points, (x,y))
    end
    return line_points
end


function k_path(p::ParamsHF, num_points; ϵ = 0)
    g1, g2 = bravais_vectors(p)
    K1 =  κ(g1, g2)
    K2 =  κ(g2, g1)
    # Initialize an array to store the line points
    # line_points = Array{Tuple{Float64, Float64}}(undef, num_points+1)
    # Compute the line points and store them in the array
    line_points = []
    for i in 0:2(√3/2+1/2) * num_points+1
        if i ≤ 1/2* num_points
        x = -K1[1]
        y = K1[2] + 2i *(1+ϵ) * K2[2]/num_points

        elseif i >1/2*num_points && i ≤ (√3/2+1/2) * num_points
        x = -K1[1] + (i-1/2*num_points) * K1[1]/(√3/2*num_points)
        y = 0.0

        elseif i >(√3/2+1/2) * num_points && i ≤ (2√3/2+1/2) * num_points
        x = 0.0 + (i-(√3/2+1/2)*num_points) * K1[1]/(√3/2*num_points)
        y = 0.0
        else
        x = K1[1]
        y = 0 + (i-(2√3/2+1/2) * num_points) * K2[2]/(num_points/2)
        # line_points[i+1] = (x, y)
        end
    push!(line_points, (x,y))
    end
    return line_points
end

# function k_path(symbols, points)

# function k_path(k1, k2, points)
#     klist = []
#     for 1:points
#        push!(klist, (k1[1] + i/num_points * (k2[1]-k1[1]), k1[2] + i/num_points * (k2[2]-k1[2])))
#     end
#     klist
# end


# Run  es = bands(101; M = 0); plotbands(es)

#-------------------------------------------------------------------------------------------
#                            Berry connection
#-------------------------------------------------------------------------------------------

# Pseudo code
# provide a list of Ks
# compute the eigenvalues at each of these Ks and order them by their eigenvalues
# store them
# for each of these Ks compute the derivative in x and y direction
# form the bracket

function berry_connection(k0; num_central_diff_points = 5, kws...)
    # note do to the form of the band structure of the heavy-fermion model in the
    # Gamma -> K path there are no band crossings, so we do not need to disentangle bands if
    # we look into this interval
    
    # generate a list of k-vectors centered at k0 along the x and y directions
    k0xs, k0ys = _k_local_paths(k0, num_central_diff_points) # we use 5 values by default
    # generate the e0s and phi0s. Diagonalization
    phis_arr_x, phis_arr_y = _wavefunct(num_central_diff_points, k0xs, k0ys; kws...)
    # compute the derivative centered tat k0 in the x and y direction of an eigenstate i
    dxun, dyun = 
        _derivative_allmodes(k0xs, k0ys, phis_arr_x, phis_arr_y, num_central_diff_points)
    # make the brackets Ax, Ay = <u_m(k0)|(\partial_x, \partial_y)| u_n(k0)>
    um = phis_arr_x[:,:, Int((num_central_diff_points-1)/2+1)] #wavefunctions at k0
    num_modes = size(phis_arr_x,1)
    # compute the connections
    Ax, Ay = _connections(num_modes, um, dxun, dyun)
    
end

#generates the k_vecs along the x and y directions around a k-point k_0
function _k_local_paths(k0, num_points; η = 6, theta =  1.09)
    if iseven(num_points)
        num_points += 1
    else nothing end 
    # Initialize an array to store the line points
    line_points_x = Array{Tuple{Float64, Float64}}(undef, num_points)
    line_points_y = Array{Tuple{Float64, Float64}}(undef, num_points)
    # steps
    K = rot_top([4π/3a0; 0], theta/2) 
    K´= rot_bot([4π/3a0; 0], -theta/2) 
    ΔK = norm(K-K´)
    # define the step used in the evaluation of the central differences scheme
    delta = ΔK/η
    # Compute the line points and store them in the array
    for i in 1:num_points
        x1 = k0[1] -(num_points-1)/2 * delta + (i - 1) * delta
        y1 = k0[2]
        line_points_x[i] = (x1, y1)
        x2 = k0[1] 
        y2 = k0[2] -(num_points-1)/2 * delta + (i - 1) * delta
        line_points_y[i] = (x2, y2)
    end
    return line_points_x, line_points_y
end

# compute the array of wavefunctions ordered in energy
function _wavefunct(num_central_diff_points, k0xs, k0ys; kws...)
    phis_arr_x = zeros(ComplexF64, 6, 6, num_central_diff_points)
    phis_arr_y = zeros(ComplexF64, 6, 6, num_central_diff_points)
    for l in 1:num_central_diff_points
        e0, phi0 = eigen(Matrix(h_hf([k0xs[l][1], k0xs[l][2]]; kws...)))
        phis_arr_x[:,:,l] = phi0
        e0, phi0 = eigen(Matrix(h_hf([k0ys[l][1], k0ys[l][2]]; kws...)))
        phis_arr_y[:,:,l] = phi0
    end
    return phis_arr_x, phis_arr_y
end

function _derivative_allmodes(k0xs, k0ys, phis_arr_x, phis_arr_y, num_central_diff_points)
    num_modes = size(phis_arr_x,1)
    dphisx = zeros(ComplexF64, num_modes, num_modes)
    dphisy = zeros(ComplexF64, num_modes, num_modes)
    for i in 1:num_modes #iterate over the modes
        dphisx[:,i,:], dphisy[:,i,:] = _derivative_singlemode(k0xs, k0ys, phis_arr_x[:,i,:], 
            phis_arr_y[:,i,:], k0, num_central_diff_points)
    end
    return dphisx, dphisy
end

function _derivative_singlemode(kvectx, kvecty, phisx, phisy, k0, num_points)
    num_modes = size(phisx,1)
    dphisx = zeros(ComplexF64, num_modes)
    dphisy = zeros(ComplexF64, num_modes)
    # Find the index of k0 in the kvect array
    k_index = Int((num_points-1)/2 + 1)
    # Compute the 2D derivatives using central difference method
    dphisx = (phisx[:,k_index+1] .- phisx[:,k_index-1]) ./
        (kvectx[k_index+1][1] - kvectx[k_index-1][1])
    dphisy = (phisy[:,k_index+1] .- phisy[:,k_index-1]) ./ 
        (kvecty[k_index+1][2] - kvecty[k_index-1][2])
    return dphisx, dphisy
end

function _connections(num_modes, um, dxun, dyun)
    Ax = zeros(ComplexF64, num_modes, num_modes)
    Ay = zeros(ComplexF64,  num_modes, num_modes)
    for i in 1:num_modes
        Ax[i,:] = um[i,:]' * dxun
        Ay[i,:] = um[i,:]' * dyun
    end
    return Ax, Ay
end