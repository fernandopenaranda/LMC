"""
σ_xxx/σ_xx in units of [1/T]

see `magneto_conductivity` and `in_plane_bindependent_conductivity`
"""

function mc_ratio(i::Symbol, j::Symbol, k::Symbol, p, of; T = 10, τ = 100, evals = 10, Ω_contr = true, omm_contr = true)
    τ *= femto_to_seconds
    M, xmin, xmax = int_boundaries(p)
    integrand(q, num) = k_mc_ratio(i, j, k, p, q, of; T = T, τ = τ, Ω_contr = Ω_contr, omm_contr = omm_contr, num = num)
    σxxx(q) = integrand(q, true)
    σxx(q) = integrand(q, false)
    
    val, err = Cubature.hcubature(σxxx, [0, xmin[2]], [xmax[1]/2, xmax[2]]; reltol=1e-5, abstol=0, maxevals = evals);
    val2, err2 = Cubature.hcubature(σxx, [0, xmin[2]], [xmax[1]/2, xmax[2]]; reltol=1e-5, abstol=0, maxevals = evals)
    return val, val2, val/val2
end

""" 
dhx and dhy do not include a momentum dependence
these have to change if we include a rz in the hamiltonian (for instance due to a sztz mass)
"""
function k_mc_ratio(i, j, k, p, q, of; T = 0, τ = 1e-15, Ω_contr = true, omm_contr = true, num = true)
    err_dirs(i,j,k)
    h = hf_valley_spin_hamiltonian(p, q, of)
    ϵs, ψs = eigen(Matrix(h))
    dhx = dhf_hamiltonian(p, :x, of)
    dhy = dhf_hamiltonian(p, :y, of)
    dhi = ifelse(i == :x, dhx, dhy)
    rzmat = rz(ψs, p, q, of) * a0 # <-! # is the a0 needed !!! units and ξ???
    
    Cup = 1e3 * 2π * τ/(4π^2*ang_to_m^2)
    Cdown = 2π * τ* ħ_mev_s/(4π^2*ang_to_m^2) 
    # 1e3 is needed to express the electron charge in meV so it cancels 1e3 Kb T in the denominator
    # ħ_mev_s comes from the B independent term
    if num ==true
        σ = Cup * k_linear_mr_integrand(i, j, k, ϵs, ψs, rzmat, dhx, dhy, p.μ, 
            T, Ω_contr = Ω_contr, omm_contr = omm_contr)
    else
        σ = Cdown * k_in_plane_bindependent_conductivity_integrand(i, ϵs, ψs, dhi, p.μ, T)        
    end
     return σ
end



##### ATLAS CALCULATIONS

function mc_map_correlated(folder, PID, job_id, a, b, c, p, mus, ofmats, ns; rs = "", kws...)
    start_time = time_ns()
    val = mc_ratio(a, b, c, ParamsHF(p, μ = mus), ofmats; kws...)
    which_current = string(a)*string(b)*string(c)*"_"*rs*"_"
    save_data_atlas_mr(which_current, folder, PID, job_id, [ns], [val],  rs)
    elapsed_time = (time_ns() - start_time) / 10^9
    println("Elapsed time: ", elapsed_time, " seconds")
    return ns, val
end

function save_data_atlas_mr(which_current, folder, PID, job_id, array1, array2, rs = "")
    df_array1 = DataFrame(omegas = vec(array1))
    df_array2 = DataFrame(omegas = vec(array2))
    # Define file names
    str = "/scratch/ferpe/HFCorrelated/cluster_temp/src" * "/output/Data" * "/" * "MR" * string(PID) * "-" *  folder
    if !isdir(str)
        mkpath(str)
    end
    which_current *= string(job_id)
    file_array1 = str * "/" * which_current * "_ns_"*rs*".csv"
    file_array2 = str * "/" * which_current * "_mr_"*rs*".csv"
    CSV.write(file_array1, df_array1)
    CSV.write(file_array2, df_array2)
    # Write to CSV files
    println("Files saved:")
    println("  - ", file_array1)
    println("  - ", file_array2)
end

function save_eval_presets(rs, evals, T, PID, job_id, folder)
    filename = "/scratch/ferpe/HFCorrelated/cluster_temp/src" * "/output/Data" * "/" * "MR" * string(PID) * "-" *  folder
    if job_id == 1
        df = DataFrame(gs = rs, evals = evals, T = T)
        CSV.write("/eval_presets" * filename, df)
    else nothing end
end 