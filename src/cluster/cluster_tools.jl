const SCRIPT_DIR = normpath(joinpath(@__DIR__)) * "/src/cluster/bash/"

function script_path(name::String)
    path = joinpath(SCRIPT_DIR, name)
    isfile(path) || error("Script not found: $name")
    return path
end

function slurm_submit_interpolations(;Ezmin = -6, Ezmax = 6, evals = 10, N = 7, eta = 0.05; 
    dryrun=false)
    script = script_path("run_interpolations.sh")
    cmd = `sbatch $script $Ezmin $Ezmax $evals $N $eta`
    dryrun && return cmd
    run(cmd)
end