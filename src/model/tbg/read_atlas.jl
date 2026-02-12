function provide_folder_get_observables_atlas_reshuffled(folder, rs = "rs")
    str = proj_folder*"/Paper/ClusterTBG/"*string(folder)
       es_path = str * "/joined_es_"*rs*".csv"
       mus_path = str * "/joined_mus_"*rs*".csv"
       ns_path = str * "/joined_ns_"*rs*".csv"
       ofmats_path = str * "/joined_ofmats_"*rs*".csv"
    #    s = return_variables_reshuffled(es_path,mus_path, ns_path, ofmats_path )
       Self_consistent_data(degen_decission(s), s[2], s[3], s[4])
       Self_consistent_data(s[1], s[2], s[3], s[4])
       
end

"noise instabilities over deg groundstates reduction"
function degen_decission(s)
    nss = [[s[1][i][fill] for i in 1:length(s[2])] for fill in 1:8];
    intnss = [interpolate_despike(s[4], nss[fill], 0.01) for fill in 1:8];
    return [getindex.(intnss, i) for i in 1:length(intnss[1])]
end

function substrateoff(PID::Array)
    s = [substrateyesornot(p) for p in PID]
    numbers[findall(x-> x == false, s)]
end

function interpolate_despike(x,y, threshold)
    good = detect_outliers(y, threshold)
    xg = x[good]
    yg = y[good]
    itp = linear_interpolation(xg, yg, extrapolation_bc = Line())
    return itp.(x)
end

function detect_outliers(data::AbstractVector{<:Real}, threshold::Real)
    n = length(data)
    good = trues(n)
    for i in 2:n-1
        avg_neighbors = (data[i-1] + data[i+1]) / 2
        if abs(data[i] - avg_neighbors) > threshold
            good[i] = false
        end
    end
    return good
end

function despike_forward!(d, threshold)
    data = zeros(length(d))
    n = length(data)
    data[1] = d[1]
    for i in 2:n-1
        left  = data[i-1]      # already filtered
        mid   = d[i]
        right = d[i+1]      # original (not yet filtered)

        avg_neighbors = (left + right) / 2
        deviation = abs(mid - avg_neighbors)

        if deviation > threshold && abs(right)-abs(mid) > threshold && abs(left)-abs(mid) > threshold || deviation > 10threshold
        data[i] = avg_neighbors
        else 
            data[i] = d[i]
        end
    end

    return data
end


substrateyesornot(PID::Number) = substrateyesornot(string(PID))
function substrateyesornot(PID::String)
    file = proj_folder*"/Paper/ClusterTBG/"*PID*"/"*"params.csv"
    if isfile(file)
        l = CSV.read(file, DataFrame)
        if :sigmaz in names(l) && :sigmazlayerz in names(l)
            return ifelse(abs(l.sigmaz[1] + l.sigmazlayerz[1]) > 1e-3, true, false)
        else 
            return ifelse(abs(l.sigmaz[1]) > 0, true, false)
        end
    else nothing end
end
readparams(PID::Number) = readparams(string(PID))
function readparams(PID::String)
    file = proj_folder * "/Paper/ClusterTBG/" * PID * "/" * "params.csv"
    if isfile(file)
        l = CSV.read(file, DataFrame)
        return l
    end
end
"""
self-consistent magnetoconductivity sweep with filling
opts: "vp", "qah"
"""
function drude_sweep(folder, phase = "vp"; kws...)
    s = provide_folder_get_observables_atlas_reshuffled(folder, phase)
    l = readparams(folder)
    pimported = ParamsHF(paramsHF(1.05, 1,  1),
        μ = 0,
        nf = l.nf[1],
        λ = l.λ[1],
        M = l.M[1],
        v = l.v[1],
        vp = l.vp[1],
        γ =l.γ[1],
        sigmaz = l.sigmaz[1],
        sigmazlayerz = l.sigmazlayerz[1],
        U1 = l.U1[1],
        U2 = l.U2[1],
        J = 0*l.J[1],
        VP =l.VP[1],
        twovalleystwospins = l.twovalleystwospins[1])
        drude_sweep(pimported, s; kws...)
end

function lmc_orb_sweep(folder, phase = "vp"; kws...)
    s = provide_folder_get_observables_atlas_reshuffled(folder, phase)
    l = readparams(folder)
    pimported = ParamsHF(paramsHF(1.05, 1,  1),
        μ = 0,
        nf = l.nf[1],
        λ = l.λ[1],
        M = l.M[1],
        v = l.v[1],
        vp = l.vp[1],
        γ =l.γ[1],
        sigmaz = l.sigmaz[1],
        sigmazlayerz = l.sigmazlayerz[1],
        U1 = l.U1[1],
        U2 = l.U2[1],
        J = 0*l.J[1],
        VP =l.VP[1],
        twovalleystwospins = l.twovalleystwospins[1])
        lmc_orb_sweep(pimported, s; kws...)
end

function lmc_spin_sweep(folder, phase = "vp"; kws...)
    s = provide_folder_get_observables_atlas_reshuffled(folder, phase)
    l = readparams(folder)
    pimported = ParamsHF(paramsHF(1.05, 1,  1),
        μ = 0,
        nf = l.nf[1],
        λ = l.λ[1],
        M = l.M[1],
        v = l.v[1],
        vp = l.vp[1],
        γ =l.γ[1],
        sigmaz = l.sigmaz[1],
        sigmazlayerz = l.sigmazlayerz[1],
        U1 = l.U1[1],
        U2 = l.U2[1],
        J = 0*l.J[1],
        VP =l.VP[1],
        twovalleystwospins = l.twovalleystwospins[1])
        lmc_spin_sweep(pimported, s; kws...)
end

function drude_sweep(p::ParamsHF, s::Self_consistent_data; kws...)
    d_arr = []
    for (it,ν) in enumerate(s.ns)
        if it % 5 == 0
            d = drude_conductivity(
                xx_drude_presets(ParamsHF(p, μ = s.mus[it]),  diagm(s.ofmats[it]); kws...))
            append!(d_arr, d)
            println("progress: ", it/length(s.ns)," |  val: ", d)
        end
    end
    return s.ns, d_arr
end

 
lmc_orb_sweep(p::ParamsHF, s::Self_consistent_data; kws...) = 
    mr_orb_sweep(linear_magneto_conductivity_orbital, p, s; kws...)

lmc_spin_sweep(p::ParamsHF, s::Self_consistent_data; kws...) = 
    mr_spin_sweep(linear_magneto_conductivity_spin, p, s; kws...)

function mr_orb_sweep(obs, p::ParamsHF, s::Self_consistent_data; kws...)
    mc_arr = []
    for (it,ν) in enumerate(s.ns)
        if it == 130
            mc = obs(
                    xxx_lmc_presets(ParamsHF(p, μ = s.mus[it]),  diagm(s.ofmats[it]); kws...))
            append!(mc_arr, mc)
            println("progress: ", it/length(s.ns)," |  val: ", mc)
        end
    end
    return s.ns, mc_arr
end


function mr_spin_sweep(obs, p::ParamsHF, s::Self_consistent_data; kws...)
    mc_arr = []
    for (it,ν) in enumerate(s.ns)
        mc = obs(
            Planar_σijk_presets_spin([1I,1I,1I], 
                xxx_lmc_presets(ParamsHF(p, μ = s.mus[it]), 
            diagm(s.ofmats[it]); kws...)))
        append!(mc_arr, mc)
        println("progress: ", it/length(s.ns)," |  val: ", mc)
    end
    return s.ns, mc_arr
end


# using CSV
# using DataFrames
# # reshapedata(PID) #run
# using DataFrames, LinearAlgebra, CSV
# using PooledArrays

#PID =  "MR2567402-2281278"
# reshapemrdata(PID::Number) = reshapedata(string(PID))
# function reshapemrdata(PID::String)
#     str0 = "/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Projects/HeavyFermion_Optics/MRDATA" * "/" * PID
#     matching_mrs = filter(x -> endswith(basename(x), "_mr_qah.csv") || endswith(basename(x), "_mr_vp.csv"), readdir(str0, join=true))
#     matching_ns = filter(x -> endswith(basename(x), "_ns_qah.csv") || endswith(basename(x), "_ns_vp.csv"), readdir(str0, join=true))
#     nlist = []
#     mrxxxlist = []
#     mrxxlist = []
#     mrratlist = []
#     for mr in matching_mrs
#         σxxx, σxx, rat = parse_pooled_tuple(CSV.read(mr, DataFrame).omegas)
#         push!(mrxxxlist, σxxx)
#         push!(mrxxlist, σxx)
#         push!(mrratlist, rat)
#     end

#     for n in matching_ns
#         push!(nlist, CSV.read(n, DataFrame).omegas[1])
#     end        
#     # save_arrays(str0 * "/joined_ns.csv", nlist)
#     # save_arrays(str0 * "/joined_mr.csv", mrlist)
#     return nlist, mrxxxlist, mrxxlist, mrratlist
# end

# readpresets(PID::Number) = readpresets(string(PID))
# function readpresets(PID::String)
#     l = CSV.read("/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Projects/HeavyFermion_Optics/MRDATA" * "/" * PID * "/eval_presets.txt", DataFrame)
#     return l
# end

# function parse_pooled_tuple(pv::PooledArrays.PooledVector{String})
#     str = strip(pv[1], ['(', ')', '"'])  # Get the string and clean it
#     parts = split(str, ",")
#     a = parse(Float64, strip(parts[1]))
#     b = parse(Float64, strip(parts[2]))
#     c = parse(Float64, strip(parts[3]))
#     return a, b, c
# end