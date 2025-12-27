# job_id = parse(Int, ARGS[1])
# jobs_num = parse(Int, ARGS[2])
# PID = ARGS[3]
# Ezmin = ARGS[4]
# Ezmax = ARGS[5]
# evals = ARGS[6]
# N = ARGS[7]
# η = ARGS[8]

# print("Starting...")

# using LMC, CSV

# Ez = range(Ezmin, Ezmax, step = (Ezmax-Ezmin)/jobs_num)[job_id+1]
# estimated_bound_width = 10

# cpt = Computation_params(estimated_bound_width, evals, η)
# p = Params_rhombohedral(1, 0, 3160, 390, -20, 315, 44, 0, 0)
# intp = Interpolated_params(N, p, [Ez], cpt)

# data_folder = pwd() * "/Data/Interpolations/" * string(PID) * "/" * string(job_id) * "/"
# save_to_csv(intp, data_folder * "presets.csv")

# Ezs, ϵ_mat, int_dos_mat, int_n_mat = interpolated_dos_ns_Ez(intp);
# @save data_folder * "interpolateddata.jld" Ezs ϵ_mat int_dos_mat int_n_mat

# str = pwd() * "/slurm-" * string(PID) * "." * string(job_id)
# str_dest = pwd() * "/Data" * "/" * string(PID) * "/" * string(job_id) * "/"
# mv(str * ".out", str_dest * "output.out")
# mv(str * ".err", str_dest * "error.err")