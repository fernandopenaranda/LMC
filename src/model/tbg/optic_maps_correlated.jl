#@time self_consistent_shift_map(:REAL, :x, :x, :x, 
#    ParamsHF(p, sigmaz = 2, vafek = 2), collect(0:.2:15), collect(-4:.1:4), η = 0.5, evals = 10000)
struct Self_consistent_data
    ofmats
    mus
    es
    ns
end

self_consistent_shift_map_correlated(folder, PID, part, a, b, c, p, sc_object::Self_consistent_data, ωlist; kws...) = 
    self_consistent_shift_map_correlated(folder, PID, part, a, b, c, p, sc_object.mus, diagm.(sc_object.ofmats), sc_object.ns, ωlist; kws...)

function self_consistent_shift_map_correlated(folder, PID, part, a, b, c, p, mus, ofmats, ns, ωlist; rs = "", gs = "", kws...)
    start_time = time_ns()
    mat1 = zeros(Float64, length(ωlist), length(ns))
    count = 1;
    for i in 1:length(ns)
        print("  |    Status: $(round(100*count/length(ns), digits = 1))%:      ")
        vals = shift_current(a, b, c, ParamsHF(p, μ = mus[i], n= ns[i]), ofmats[i], ωlist, part; kws...)[2]
        mat1[:,count] .= vals
        count +=1
    end
        which_current = "shift"
    if part == :REAL
        which_current *= "_real_"
    else 
        which_current *= "_imag_"
    end
	which_current *= string(a)*string(b)*string(c)*rs
    save_data_atlas(which_current,folder, PID, ns, ωlist, mat1,  gs)
    elapsed_time = (time_ns() - start_time) / 10^9
    println("Elapsed time: ", elapsed_time, " seconds")
    return ns, ωlist, mat1
end

self_consistent_injection_map_correlated(folder, PID, part, a, b, c, p, sc_object::Self_consistent_data, ωlist; kws...) = 
    self_consistent_injection_map_correlated(folder, PID, part, a, b, c, p, sc_object.mus, diagm.(sc_object.ofmats), sc_object.ns, ωlist; kws...)

function self_consistent_injection_map_correlated(folder, PID, part, a, b, c, p, mus, ofmats, ns, ωlist;rs = "",gs ="",  kws...)
    start_time = time_ns()
    mat1 = zeros(Float64, length(ωlist), length(ns))
    count = 1;
    for i in 1:length(ns)
        print("  |    Status: $(round(100*count/length(ns), digits = 1))%:      ")
        vals = injection_current(a, b, c, ParamsHF(p, μ = mus[i], n= ns[i]), ofmats[i], ωlist, part; kws...)[2]
        mat1[:,count] .= vals
        count +=1
    end
        which_current = "inj"
    if part == :REAL
        which_current *= "_real_"
    else 
        which_current *= "_imag_"
    end
	which_current *= string(a)*string(b)*string(c)*rs
    save_data_atlas(which_current,folder, PID, ns, ωlist, mat1, gs)
    elapsed_time = (time_ns() - start_time) / 10^9
    println("Elapsed time: ", elapsed_time, " seconds")
    return ns, ωlist, mat1
end

provide_folder_get_observables(folder::Number) = 
    provide_folder_get_observables(string(folder))

function provide_folder_get_observables(folder)
    str0 = "/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Projects/HeavyFermion_Optics/src/data_atlas_selfconsistency/"
    str = str0 * folder
    es_path = str * "/joined_es.csv"
    mus_path = str * "/joined_mus.csv"
    ns_path = str * "/joined_ns.csv"
    ofmats_path = str * "/joined_ofmats.csv"
    s = return_variables(es_path,mus_path, ns_path, ofmats_path )
    Self_consistent_data(s[1], s[2], s[3], s[4])
end

function provide_folder_get_observables_atlas(folder)
    str = pwd() * "/Data" * "/" * string(folder)
    es_path = str * "/joined_es.csv"
    mus_path = str * "/joined_mus.csv"
    ns_path = str * "/joined_ns.csv"
    ofmats_path = str * "/joined_ofmats.csv"
    s = return_variables(es_path,mus_path, ns_path, ofmats_path )
    Self_consistent_data(s[1], s[2], s[3], s[4])
end

# function provide_folder_get_observables_atlas_reshuffled(folder, rs = "rs")
#  str = pwd() * "/Data" * "/" * string(folder)
#     es_path = str * "/joined_es_"*rs*".csv"
#     mus_path = str * "/joined_mus_"*rs*".csv"
#     ns_path = str * "/joined_ns_"*rs*".csv"
#     ofmats_path = str * "/joined_ofmats_"*rs*".csv"
#     s = return_variables_reshuffled(es_path,mus_path, ns_path, ofmats_path )
#     Self_consistent_data(s[1], s[2], s[3], s[4])
# end

function return_variables(es_path,mus_path, ns_path, ofmats_path )
    ns_unordered = flatten(parse_str_to_arr.(import_array(ns_path)))
    indices = sortperm(ns_unordered)
    ns = ns_unordered[indices]
    es = import_array(es_path)[:][indices]
    mus = import_array(mus_path)[indices]
    ofmats = import_array(ofmats_path)
    ofs = [parse_str_to_arr(of[3:end-2]) for of in ofmats]
    return ofs[indices], mus, es, ns
end


function return_variables_reshuffled(es_path,mus_path, ns_path, ofmats_path )
    ns_unordered = import_array(ns_path)
    indices = sortperm(ns_unordered[:])
    ns = ns_unordered[indices]
    es = import_array(es_path)[:][indices]
    mus = import_array(mus_path)[:][indices]
    ofmats = import_array(ofmats_path)
    ofs = [parse_str_to_arr(of[3:end-2]) for of in ofmats]
    return ofs[indices], mus, es, ns
end


index_to_closestval(arr, val) = argmin(abs.(arr .- val))
parse_str_to_arr(str) = 
           map(x -> parse(Float64, x), split(strip(str, ['[', ']']), ","))

four_op_vsfilling(p, sc_obj::Self_consistent_data; kws...) = 
    four_op_vsfilling(p, sc_obj.mus, diagm.(sc_obj.ofmats); kws...)


function save_data_atlas(which_current, folder, PID, array1, array2, matrix1, rs = "")
    df_array1 = DataFrame(omegas = vec(array1))
    df_array2 = DataFrame(omegas = vec(array2))

    # Define file names
    str = "/scratch/ferpe/HFCorrelated/cluster_temp/src" * "/output/Data" * "/" * string(folder) * "/" * string(PID) * "/"
    file_array1 = str * which_current * "array1"*rs*".csv"
    file_array2 = str * which_current * "array2"*rs*".csv"
    file_matrix1 = str * which_current * "mat"*rs*".csv"
    CSV.write(file_array1, df_array1)
    CSV.write(file_array2, df_array2)
    CSV.write(file_matrix1,  Tables.table(matrix1), writeheader=false)

    # Write to CSV files
    println("Files saved:")
    println("  - ", file_array1)
    println("  - ", file_array2)
    println("  - ", file_matrix1)
end


import_array(file_arrays) =  Matrix(CSV.read(file_arrays, DataFrame))
