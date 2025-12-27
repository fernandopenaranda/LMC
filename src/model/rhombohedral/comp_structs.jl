@with_kw struct Computation_params
    estimated_bound_width::Number
    evals::Number
    η::Float64
end

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

function save_to_csv(ip::Interpolated_params, filename::String)
    df = DataFrame([to_dict(ip)])
    CSV.write(filename, df)
end