@with_kw struct Computation_params
    estimated_bound_width::Number
    evals::Number
    η::Float64
    λ::Number = 1e5## penalty for not conserved number of particles
    iterations::Int = 10
    random_guesses::Int = 20
end

Computation_params(estimated_bound_width, evals, η) = 
    Computation_params(estimated_bound_width, evals, η, 1e5, 10,20)

@with_kw struct Interpolated_params
    N::Int # Number of layers
    p::Params_rhombohedral
    Ezlist::Array # Displacement fields
    cpt::Computation_params
end


function to_dict(cpt::Computation_params)
    return Dict(
        :estimated_bound_width => cpt.estimated_bound_width,
        :evals => cpt.evals,
        :η => cpt.η)
end

function to_dict(ip::Interpolated_params)
    d = Dict{Symbol,Any}()
    d[:N] = ip.N
    # Params_rhombohedral (store as readable string)
    d[:p] = repr(ip.p)
    # Ezlist as comma-separated string
    d[:Ezlist] = join(ip.Ezlist, ",")
    # flatten Computation_params
    for (k, v) in to_dict(ip.cpt)
        d[k] = v
    end
    return d
end

# function save_to_csv(ip::Union{Interpolated_params, Phase_diagram_params}, filename::String)
#     df = DataFrame([to_dict(ip)])
#     CSV.write(filename, df)
# end

@with_kw struct Drude_params
    N::Int # Number of layers
    p::Params_rhombohedral
    Ezlist::Array # Displacement fields
    cpt::Computation_params
end

@with_kw struct Phase_diagram_params
    N::Int # Number of layers
    p::Params_rhombohedral
    Ezsteps::Number # Displacement fields
    νmin::Float64
    νmax::Float64
    νpoints::Int
    int_model::Symbol
    U::Number
    J::Number
    cpt::Computation_params
end

@with_kw struct Observable_computation
    job_id
    jobs_num
    PID
    interpPID  # PID of the process that generated the interpolated data
    phasediagPID    # PID of the process that generated the hartree mus and nss out of interpolated data with PID interpPID
    evals
    T # temperature a <= eta/kb ~ 1K, since this is really a T=0 limit and in the Phase diagram T = 0 and there are only broadening effects
    tau
    which_observable
end