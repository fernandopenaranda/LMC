for (i, a) in enumerate(ARGS)
    @info "ARGS[$i] = $a"
end

job_id = parse(Int, ARGS[1])
jobs_num = parse(Int, ARGS[2])
PID = ARGS[3]
dataPID = ARGS[4]     # PID of the process that generated the interpolated data
Ezsteps = parse(Int, ARGS[5])
nu_min = parse(Float64, ARGS[6])
nu_max = parse(Float64, ARGS[7])
nu_points = parse(Int, ARGS[8])
U = parse(Float64, ARGS[9])
J = parse(Float64, ARGS[10])
λ =parse(Float64, ARGS[11])
η = parse(Float64, ARGS[12])
estimated_bound_width = parse(Int, ARGS[13])
iterations = parse(Int, ARGS[14])
int_model = Symbol(ARGS[15])
random_guesses = parse(Int, ARGS[16])
evals = parse(Int, ARGS[17])

println("Importing interpolating data")

using Pkg
println(Pkg.status())

using JLD2, CSV, DataFrames, Interpolations, LMC
filestring = pwd() * "/Data/Interpolations/" * string(dataPID) * "/" * 
    string(job_id) 
@load filestring * "/interpolateddata.jld" Ezs ϵ_mat int_dos_mat int_n_mat
@load filestring * "/presets.csv" intp

#ensure that there is always a root in findroot see Eμαs findzero
hard_up_bound = maximum(int_n_mat[1].itp)
hard_low_bound = minimum(int_n_mat[1].itp)

if nu_min < hard_low_bound  
    nu_min = hard_low_bound
else nothing end
if nu_max> hard_up_bound
    nu_max = hard_up_bound
else nothing end


print("Building structs...")
cpt = Computation_params(estimated_bound_width, evals, η, λ, iterations,
    random_guesses)

PD_presets = Phase_diagram_params(intp.N, intp.p, Ezsteps, nu_min, nu_max, nu_points,
    int_model, U, J, cpt)

print("Computation...")

@time μs, ns = Emin_nαs(int_dos_mat[1], int_n_mat[1], PD_presets);

print("Storing...")

new_data_folder = pwd() * "/Data/PhaseDiagrams/" * string(PID) * "/" * string(job_id)
mkpath(new_data_folder)

@save new_data_folder * "/presets.csv" PD_presets
write(new_data_folder * "/output.txt", "Using interpolated data from $(dataPID)")

nu_list = collect(range(nu_min, nu_max, step = (nu_max-nu_min)/nu_points))
mus = μs[1][1]
nss = ns[1][1]
@save new_data_folder * "/pd_data.jld" nu_list mus nss Ezs

print("Success!")

str = pwd() * "/slurm-" * string(PID) * "." * string(job_id)
mv(str * ".out", new_data_folder * "/output.out")
mv(str * ".err", new_data_folder * "/error.err")