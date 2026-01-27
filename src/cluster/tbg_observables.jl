#=
data generated in our previous paper based on a self-consistent hartree fock solution
of the topological heavy fermion model on the interacting bandstructures with SU(4)
are passed to the optical conductivity functions in Optics_in_the_length_gauge. LMC
is a wrapper.
Data generated from the THFM is stored in Paper/ClusterTBG
    observable: Drude, LMC_orbital, LMC_spin or QAH
    phases = "vp", "qah"
=#
for (i, a) in enumerate(ARGS)
    @info "ARGS[$i] = $a"
end

job_id = parse(Int, ARGS[1])
jobs_num = parse(Int, ARGS[2])
PID = ARGS[3]
folder = ARGS[4]
phase = ARGS[5]
evals = parse(Int, ARGS[6])
T = parse(Float64, ARGS[7]) # temperature a <= eta/kb ~ 1K, since this is really a T=0 limit and in the Phase diagram T = 0 and there are only broadening effects
tau = parse(Float64, ARGS[8])
which_observable = ARGS[9]  # Drude, LMC_orbital, LMC_spin, or QAH

println("Importing interpolating data")

using JLD2, CSV, DataFrames, Interpolations, LMC, Optics_in_the_length_gauge, LinearAlgebra

# import data from the hartree calculation and presets
s = LMC.provide_folder_get_observables_atlas_reshuffled(folder, phase)
l = LMC.readparams(folder)
p = ParamsHF(paramsHF(1.05, 1,  1),
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


#define the observable and its presets
if which_observable == "Drude"
    obs = drude_conductivity
    press = xx_drude_presets(ParamsHF(p, μ = s.mus[job_id+1]), 
        diagm(s.ofmats[job_id+1]), T = T, τ = tau, evals = evals)


elseif which_observable == "LMC_orbital"
    obs = linear_magneto_conductivity_orbital
    press = xxx_lmc_presets(ParamsHF(p, μ = s.mus[job_id+1]), 
        diagm(s.ofmats[job_id+1]), T = T, τ = tau, evals = evals, ϵ = 1e-7,  
        berry_contribution = true, omm_contribution = true, 
        fermi_surface = false, with_shift= false)
   
elseif which_observable == "LMC_spin"
    obs = linear_magneto_conductivity_spin
    # press = xxx_lmc_spin_presets(p::ParamsHF, of; T = 1, τ = 200, evals = 100, ϵ = 1e-7,  
    # berry_contribution = true, omm_contribution = true, fermi_surface = false,
    # with_shift= false)
    throw(ArgumentError("observable not allowed yet"))
    press = 0
   
elseif which_observable == "QAH" 
    obs = σij_anomalous_hall
    throw(ArgumentError("observable not allowed yet"))
    press = 0
    
else throw(ArgumentError("observable not allowed"))
end

print("Building structs...")

obs_val = obs(press)

print("Computation...")

print("Storing...")

new_data_folder = pwd() * "/Data/tbg_$(which_observable)/" * string(PID) * "/" * string(job_id)
mkpath(new_data_folder)

comp_struct = Observable_computation(job_id, jobs_num, PID, dataPID, 
    phasediagPID, evals, T, tau, which_observable)

n = s.ns[job_id+1]
mu = s.mus[job_id+1]
@save new_data_folder * "/presets.jld" comp_struct
@save new_data_folder * "/data.jld" n mu p

print("Success!")

str = pwd() * "/slurm-" * string(PID) * "." * string(job_id)
mv(str * ".out", new_data_folder * "/output.out")
mv(str * ".err", new_data_folder * "/error.err")