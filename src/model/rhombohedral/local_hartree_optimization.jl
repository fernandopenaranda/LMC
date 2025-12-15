#_________________________________________________________________________________________
# Local Hartree interaction
#_________________________________________________________________________________________
using Optim

"""
by introducing such SU2 Hunds interaction we have broken the arbitrarity between indices
of the four flavours. Now 1 → K↑, 2 → K'↑, 3 → K↓, 4 → K'↓.
This function computes the order parameters associated to the Hartree phase.
    σ are the spin τ valley d.o.f
"""
function rh_order_parameter(nαs)
    σ0τ0 = Δw(nαs, [1,1,1,1])/norm(nαs)
    σzτ0 = Δw(nαs, [1,1,-1,-1])/norm(nαs)
    σ0τz = Δw(nαs, [1,-1,1,-1])/norm(nαs)
    σzτz = Δw(nαs, [1,-1,-1,1])/norm(nαs)
    return σ0τ0, σzτ0, σ0τz, σzτz
end

function Δw(n, w)
    return dot(w, n)
end



""" 
`Emin_μαs(p::Planar_σijk_presets, μ0s::Vector{Int64}, μs::Union{Array,AbstractRange};
    μ = 0, U = 0, λ = 0, evals = 1e2, η = 0.05, 
    estimated_bound_width = 10, iterations = 100) `  
performs several random seed calculations of the extrema values of the 
grand potential functional and selects the set of μαs with lower energy.
Crucially the dos and the occupation are independent of μ and thus are interpolated
before calling to the optimizer (highly efficient).
The system in degenerate in spin but not in valley, however the total energy and filling is
degenerate in valley too.
The type of potential SU4 or SU2 (hunds coupling) is selected by int_model
"""

function Emin_nαs(int_dos::ScaledInterpolation, n::ScaledInterpolation, μarray; kws...)
    nαs = []
    for μ in μarray
        push!(nαs, n.(Emin_μαs(int_dos, n, μ; kws...))) #sweep in the chemical potential for a constant value
    end
    return nαs
end

function Emin_μαs(p::Planar_σijk_presets, μ::Number; 
    evals = 1e2, η = 0.05, estimated_bound_width = 20, kws...) 
    ϵ_range, int_dos = interpolated_dos(p, 2estimated_bound_width, evals = evals, η = η)
    n = interpolated_n(int_dos, ϵ_range)
    Emin_μαs(int_dos, n, μ; evals = evals, η = η, 
        estimated_bound_width = estimated_bound_width, kws...)
end


function Emin_μαs(int_dos::ScaledInterpolation, n::ScaledInterpolation, μ::Number; 
        random_guesses = 2, U = 0, J = 0, λ = 0, evals = 1e2, η = 0.05, 
        estimated_bound_width = 10, iterations = 100, int_model = :SU4) 
    
    U *= (√3/2)* (a0)^2 # over the BZ area
    J *= (√3/2)* (a0)^2 # over the BZ area
    count = 0
    μα_mat = []
    Es = []
    while count < random_guesses
        μ0s = rand(4) * 2μ
        μαs = opt_μs(μ0s, n; μ = μ, U = U, J = J, λ = λ, 
        evals = eval, η = η, estimated_bound_width = estimated_bound_width, 
            iterations = iterations, int_model = int_model)
        E = total_energy(μαs, int_dos, n, U, J, int_model) # compute the total energy
        push!(μα_mat, μαs)
        append!(Es, E)
        count +=1 
    end
    ind = findmin(Es)[2]
    # println("E: ", Es[ind])
    # println("constraint: ∑n_α - 4n0 = ", sum([n(μi)-n(μ) for μi in μα_mat[ind]]))
    # println("minimal μαs: ", μα_mat[ind])
    return μα_mat[ind]
end


function Emin_μαs(p::Planar_σijk_presets, μ0s::Array, μ::Number;
        U = 0, J = 0, λ = 0, evals = 1e2, η = 0.05, 
        estimated_bound_width = 10, iterations = 100, int_model = :SU4) 
    U *= (√3/2)* (a0)^2 # over the BZ area
    J *= (√3/2)* (a0)^2 # over the BZ area
    ϵ_range, int_dos = interpolated_dos(p; evals = evals, η = η)
    n = interpolated_n(int_dos, ϵ_range)
    μαs = opt_μs(μ0s, n; μ = μ, U = U, J = J, λ = λ, 
        evals = eval, η = η, estimated_bound_width = estimated_bound_width, 
        iterations = iterations, int_model = int_model)
    E = total_energy(μαs, int_dos, n, U, J, int_model) # compute the total energy
    println("E: ", E)
    return μαs
end

""" 
`opt_μs(p::Planar_σijk_presets, μ0s::Vector{Int64}; evals = 10, η = 0.05, kws...) = 
    opt_μs(μ0s, interpolated_dos(p; evals = evals, η = η); kws...)`
or 
`opt_μs(μ0s, es_and_int_dos; kws...) = 
    opt_μs(μ0s, interpolated_n(es_and_int_dos[2], es_and_int_dos[1]); kws...)`
or 
`opt_μs(μ0s::Vector{Float64}, n::ScaledInterpolation; μ = 0, U = 0, λ = 0,
    evals = 1e2, η = 0.05, estimated_bound_width = 10, iterations = 100)`

method of bandstructures with a Hartree local potential
here we find the extrema of the grand functional: 
                Φ/A = ∑_iE_α(μ_α) + V_Hartree - TS - μ ∑_α n_α.
Its minima in terms of μ_α a 4-dimensional array corresponding to the chemical potentials 
of each of the 4 flavors

So we take the derivatives of the grand functional which lead (even at finite T) to the 
following extrema equations F(i) = μ_i - μ + U A * ∑_β≠α n_β = 0 with i and j in 1:4. 
And reformulate them as an optimization problem so ∑_i(F(i)^2) = 0. This is the quantity
to optimize in terms of μ_i, not that n_β = ∫₀^β ρ(ϵ) dϵ. Since ρ(ϵ) the dos is the same
for the four flavours in this rigid band (local (q independent) Hartree perturbation) we
first create this function and interpolate it. We also interpolate the density of states 
for convenience μ0s are the seeds

Attention there is a field inside the Planar_σijk_presets.computation that is evals, we do
not use it here but you should change the kws to be consistent in future calculations. Also
allow for the Dos struct.
"""
opt_μs(p::Planar_σijk_presets, μ0s::Vector{Int64}; evals = 10, η = 0.05, kws...) = 
    opt_μs(μ0s, interpolated_dos(p; evals = evals, η = η); kws...)

opt_μs(μ0s, es_and_int_dos; kws...) = 
    opt_μs(μ0s, interpolated_n(es_and_int_dos[2], es_and_int_dos[1]); kws...)

function opt_μs(μ0s::Vector{Float64}, n::ScaledInterpolation;
    μ = 0, U = 0, J = 0, λ = 0, evals = 1e2, η = 0.05, estimated_bound_width = 10, 
        iterations = 100, int_model = :SU4)
    lower = -estimated_bound_width .* ones(4)
    upper = estimated_bound_width .* ones(4)
    n0 = n(μ) # number of electrons at μ, no interaction
    if int_model == :SU4
        objective = objective_su4
    else # :SU2
        objective = objective_su2
    end 
    return opt_μs(objective, n0, n, μ0s, μ, U, J, λ, lower, upper, iterations)
end
# boxed interaction limited to small μvalues determined by lower and upper and within
# the interpolation 'mulist range
function opt_μs(objective, n0, n, μ0s, μ, U, J, λ, lower, upper, iterations)
    result = optimize(μs -> objective(n, μs, μ, U, J, λ),
            lower, # lower bounds of μαs    
            upper, # upper bounds of μαs
            μ0s,   # seed
            Fminbox(BFGS()), 
            Optim.Options(g_tol = 1e-12,
                iterations = iterations,    
                store_trace = false,
                show_trace = false,
                show_every = 5,
                show_warnings = false))

    μαs = Optim.minimizer(result)
    return μαs
end

Integ(f, μαs::Array) = [Integ(f, μα) for μα in μαs]
Integ(f, μα::Number) = quadgk(f, 0, μα)[1]

""" reformulation of the root-finder problem intro an optimization 
least squares one introducing a su4 preserving MF decoupling corresponding 
to Vint = UA/2 ∑_α≠β n_α n_β"""
function objective_su4(n, μs, μ, U, J, λ)
    n0 = n(μ)
    s = 0
    inds_not_alpha(i) = findall(j-> j≠i, 1:4)
    for (i,μi) in enumerate(μs)
        s += (μi - μ +  U * sum([n(μi) for μi in μs[inds_not_alpha(i)]]))^2  
    end
    obj = s + penalty_fixed_filling(n0, n, μs, λ = λ)
    return obj
end

""" reformulation of the root-finder problem intro an optimization least squares 
one introducing a su2 preserving MF decoupling corresponding to 
    Vint = UA/2 ∑_α≠β n_α n_β +  JAu.c.(n1 − n3)(n2 − n4)
Hund coupling where n1 and n3 are opposite spins same valley"""
function objective_su2(n, μs, μ, U, J, λ)
    n0 = n(μ)
    s = 0
    inds_not_alpha(i) = findall(j-> j≠i, 1:4)
    for (i,μi) in enumerate(μs)
        s += ( μi - μ +  
               U * sum([n(μi) for μi in μs[inds_not_alpha(i)]]) +

               J * hund_coupling(i, n, μs))^2  
    end
    obj = s + penalty_fixed_filling(n0, n, μs, λ = λ)
end

function hund_coupling(i, n, μs)
    pref = ifelse(i == 1 || i == 2, 1, -1) 
    if i == 1
        return n(μs[2]) - n(μs[4])
    elseif i == 2
        return  n(μs[1]) - n(μs[3])
    elseif i == 3
        return -n(μs[2]) + n(μs[4])
    elseif i == 4
        return -n(μs[1]) + n(μs[3])
    end

end

""" the chemical potentials resulting from the optimization should obey 
that ∑n_α(U, J)  = 4n0 with n0 the non-interacting density (equal to the 
four flavours). we do this by a lagrange multiplier"""
function penalty_fixed_filling(n0, n, μs; λ = 0)
    difn = sum([n(μi)-n0 for μi in μs])
    # println("n0: ", n0, "δn: ", difn)
    return λ * difn^2
end

"""
    `interpolated_dos(p::Planar_σijk_presets; μlist = 0.0:0.1:15.0, η = 0.05, evals = 100)`
interpolated density of statse per unit of area for an arbitrary flavour. Returns a function 
that can be evaluated within bounds of μlist. ! this code paths is at 0 temperature
"""

interpolated_dos(p::Planar_σijk_presets, bounds::Number; kws...) =
    interpolated_dos(p; μlist = -bounds:bounds/100:bounds, kws...) 

function interpolated_dos(p::Planar_σijk_presets; 
        μlist = -10.0:0.1:10.0, η = 0.05, evals = 100) 
    # println("Interpolating dos...")
    A = (1e-10)^2
    int_dos(ε) = A * c_dos(p, ε, η = η, evals = evals)[2][1]
    dos = int_dos.(μlist)
    itp = interpolate(dos, BSpline(Linear()))
    # print(" Success")
    return (μlist, Interpolations.scale(itp, μlist))
end

interpolated_n(int_dos, x::Array) = 
    interpolated_n(int_dos, range(first(x), last(x); length=length(x)))
function interpolated_n(int_dos, μlist::AbstractRange)
    n(ϵ) = Integ(int_dos, ϵ)
    ns = n.(μlist)
    itp = interpolate(ns, BSpline(Linear()))
    return Interpolations.scale(itp, μlist)
end



#_________________________________________________________________________________________
# Minimize the grand functional. Only T = 0 codepath XCommmented T ≠ 0
#_________________________________________________________________________________________
"""
the grand functional per unit of area Φ/A = sum_α [ E(μ_α)] + Vint - 
evaluated at the extrema found opt_μαs.
The entropy part which does not enter in the extrema conditions is only relevant at finite T=0
T≠0 methods not implemented at the moment.
"""
function total_energy(μαs, dos, n, U, J, int_model)
    s = 0.0
    for (i,μα) in enumerate(μαs)  
       s += E(dos, μα) - μ * n(μα)
    end
    if int_model == :SU4
        potential = v_su4
    else 
        potential = v_su2 
    end # :SU2
    s += potential(U, J, n, μαs) #- T*entropy(μαs)
end

""" total energy over the occupied bands of a single flavour at a given temperature
∫_0^μα ϵ ρ(ϵ) dϵ at T = 0 or ∫_0^∞ f(ϵ, T) ϵ ρ(ϵ) dϵ at finite T
It is # momenta independent as it should.
It uses the interpolated dos as input
Only works at zero T generalized if required
"""
function E(dos, μα)
    integrand(ϵ) = ϵ * dos(ϵ)
    return Integ(integrand, μα)    
end

"""
V_su4 = U*A/2 ∑_α≠β n_α n_β
"""
function v_su4(U, J, n, μαs)
    ns = n.(μαs) 
    s  = sum(ns)
    s2 = sum(abs2, ns)   
    return  U/2 *(s^2 - s2) # this is an identity
end

function v_su2(U, J, n, μαs)
    return v_su4(U, J, n, μαs) + J * (n(μαs[1])-n(μαs[3])) * (n(μαs[2])-n(μαs[4]))
end