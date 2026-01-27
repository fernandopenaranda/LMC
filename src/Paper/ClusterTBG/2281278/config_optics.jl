job_id = parse(Int, ARGS[1])
folder = string(ARGS[2])
PID = string(ARGS[3])


job_id += 1
include("HF_Optics.jl")

# Calculation PARAMS:
evals = 10000# 10000
omegas = collect(0:0.15:30) #0.1

p = read_struct_data(folder)
s = provide_folder_get_observables_atlas(folder)

self_consistent_shift_map_correlated(folder,job_id-1 ,:REAL, :x, :x, :x, p,  [s.mus[job_id]], [diagm.(s.ofmats)[job_id]], [s.ns[job_id]], omegas,
    η = 0.5, evals = evals)

self_consistent_shift_map_correlated(folder, job_id-1,  :REAL, :y, :y, :y, p, [s.mus[job_id]], [diagm.(s.ofmats)[job_id]], [s.ns[job_id]], omegas,
    η = 0.5, evals = evals)

self_consistent_injection_map_correlated(folder, job_id-1, :REAL, :x, :x, :x, p,  [s.mus[job_id]], [diagm.(s.ofmats)[job_id]], [s.ns[job_id]], omegas,
    η = 0.5, evals = evals)

self_consistent_injection_map_correlated(folder, job_id-1, :REAL, :y, :y, :y, p,  [s.mus[job_id]], [diagm.(s.ofmats)[job_id]], [s.ns[job_id]], omegas,
    η = 0.5, evals = evals)

job_id -= 1
str = pwd() * "/slurm-" * string(PID) * "." * string(job_id)
str_dest = pwd() * "/Data" * "/" * string(folder) * "/" * string(job_id) * "/"
mv(str * ".out", str_dest * "output_opt.out",force=true)
mv(str * ".err", str_dest * "error_opt.err", force=true)
cp("/scratch/ferpe/HFCorrelated/cluster_temp/src/atlas_optics.jl", pwd() * "/Data" * "/" * string(folder)* "/" * "config_optics.jl")
