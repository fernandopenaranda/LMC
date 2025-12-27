# const SCRIPT_DIR = normpath(joinpath(@__DIR__)) * "/src/cluster/bash/"

function script_path(name::String)
    path = joinpath(proj_folder * "/cluster/bash/", name)
    isfile(path) || error("Script not found: $name")
    return path
end

function slurm_submit_interpolations(;Ezmin = -6, Ezmax = 6, 
    evals = 10, N = 7, eta = 0.05,  dryrun=false)
    lmcfolder = dirname(pathof(LMC)) * "/cluster/interpolations.jl"
    script = script_path("run_interpolations.sh")
    println("bash_file_name: ", script)
    cmd = `sbatch $script $lmcfolder $Ezmin $Ezmax $evals $N $eta`
    dryrun && return cmd
    run(cmd)
end