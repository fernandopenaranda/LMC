job_id = parse(Int, ARGS[1])
folder = string(ARGS[2])
PID = string(ARGS[3])
rs = string(ARGS[4]) # vp or qah
evals = parse(Int, ARGS[5])
T = parse(Float64, ARGS[6])

job_id += 1
include("HF_Optics.jl")
include("degen_decision.jl") # to load the qah and vp structs
include("optic_maps_correlated.jl") # import the qah and vp structs
include("mr_module.jl")

p = read_struct_data(folder)
s = provide_folder_get_observables_atlas_reshuffled(folder, rs)


mc_map_correlated(folder,PID, job_id, :x,:x,:x, p, s.mus[job_id], diagm(s.ofmats[job_id]), s.ns[job_id],
    T = T, τ = 200, evals = evals, Ω_contr = true, omm_contr = true, rs = rs)


save_eval_presets(rs, evals, T, PID, job_id, folder) # save comp presets

job_id -= 1
str = pwd() * "/slurm-" * string(PID) * "." * string(job_id)
str_dest = pwd() * "/Data" * "/" * string(folder) * "/" * string(job_id) * "/"
mv(str * ".out", str_dest * "mr_output_opt.out",force=true)
mv(str * ".err", str_dest * "mr_error_opt.err", force=true)


