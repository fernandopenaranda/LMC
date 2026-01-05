# const SCRIPT_DIR = normpath(joinpath(@__DIR__)) * "/src/cluster/bash/"

function script_path(name::String)
    path = joinpath(proj_folder * "/cluster/bash/", name)
    isfile(path) || error("Script not found: $name")
    return path
end

function slurm_submit_interpolations(;Ezmin = -6, Ezmax = 6, 
    evals = 10, N = 7, eta = 0.05, phs = 0, dryrun=false)
    lmcfolder = dirname(pathof(LMC)) * "/cluster/interpolations.jl"
    script = script_path("run_interpolations.sh")
    println("bash_file_name: ", script)
    cmd = `sbatch --wait $script $lmcfolder $Ezmin $Ezmax $evals $N $eta $phs`
    dryrun && return cmd
    run(cmd)
    # Capture stdout
    output = read(cmd, String)
    # Example output: "Submitted batch job 3228119\n"
    jobid = parse(Int, split(output)[end])
    return jobid
end
 
""" computes the mus and alphas for the 4 spin/valley flavours
selecting the energetically favored ground state in the presense
of a local Hartree interaction of the form `int_model =:SU2 or :SU4`
`interpolationsPID::String`. Remarkably we do not need to run interpolations
of the spectral densities for each valley separately for each point in J and U 
space since the interaction is local and thus the bandstructure of the 
independent valleys is unaffected by it except from a rigid shift encoded
in mus and alphas computed in this script
"""
# slurm_submit_phasediagrams(interpolationsPID::Number; kws...) = 
#     slurm_submit_phasediagrams(string(interpolationsPID); kws...)

function slurm_submit_phasediagrams(interpolationsPID::Union{String,Number};
    dryrun=false, evals = 10, nu_min=-0.25, nu_max=0.25, nu_points=10, U=10,J=-4,lambda=1e5, eta=0.05, 
    estimated_bound_width=10, iterations=10, int_mode=:SU2, random_guesses=20)
    Ezsteps = 1 # cannot change it. Same number of jobs as in interpolations
    lmcfolder = dirname(pathof(LMC)) * "/cluster/phase_diagrams.jl"
    script = script_path("run_phasediagrams.sh")
    cmd = `sbatch --wait $script $lmcfolder $interpolationsPID $Ezsteps $nu_min $nu_max $nu_points $U $J $lambda $eta $estimated_bound_width $iterations $int_mode $random_guesses $evals`
    dryrun && return cmd
    run(cmd)

    # Capture stdout
    output = read(cmd, String)
    # Example output: "Submitted batch job 3228119\n"
    jobid = parse(Int, split(output)[end])
    return jobid
end

"""using the results of the Hartree algorithm it computes a given observable 
in Drude, LMC_orbital, LMC_spin, QAH"""
function slurm_submit_observable(phasediagramPID::Union{String,Number};
    dryrun=false, evals = 10, T = 1, tau = 200, which_observable = "Drude")
    lmcfolder = dirname(pathof(LMC)) * "/cluster/observable.jl"
    script = script_path("run_observable.sh")
    cmd = `sbatch --wait $script $lmcfolder $phasediagramPID $evals $T $tau $which_observable`
    dryrun && return cmd
    run(cmd)
    # Capture stdout
    output = read(cmd, String)
    # Example output: "Submitted batch job 3228119\n"
    jobid = parse(Int, split(output)[end])
    return jobid
end