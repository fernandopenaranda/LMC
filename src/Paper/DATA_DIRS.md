#Figure 8

# pdPID = 3238850
# drudePID = 3285558
# lmcPID = 3286797
# ahePID = 3291639
# spinPID = 3241713#3251898 #3241713


# pdpath = "/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Packages/LMC/src/Paper/ClusterDATA/LMC/$(pdPID)_merged_data.jld"
# pdpresetpath = "/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Packages/LMC/src/Paper/ClusterDATA/LMC/$(pdPID)_merged_presets.jld"
# drudepath = "/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Packages/LMC/src/Paper/ClusterDATA/LMC/$(drudePID)_merged_data.jld"
# drudepresetpath = "/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Packages/LMC/src/Paper/ClusterDATA/LMC/$(drudePID)_merged_presets.jld"
# lmcpath = "/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Packages/LMC/src/Paper/ClusterDATA/LMC/$(lmcPID)_merged_data.jld"
# lmcpresetpath = "/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Packages/LMC/src/Paper/ClusterDATA/LMC/$(lmcPID)_merged_presets.jld"
# spinpath = "/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Packages/LMC/src/Paper/ClusterDATA/LMC/$(spinPID)_merged_data.jld"
# spinpresetpath = "/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Packages/LMC/src/Paper/ClusterDATA/LMC/$(spinPID)_merged_presets.jld"


# ahepath = "/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Packages/LMC/src/Paper/ClusterDATA/LMC/$(ahePID)_merged_data.jld"
# ahepresetpath = "/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Packages/LMC/src/Paper/ClusterDATA/LMC/$(ahePID)_merged_presets.jld"


# Cluster code to genererate data  of Fig 8
# using LMC
# include("julia_exec.jl")
# pids_int = slurm_4_interpolations()
# pids_pd = slurm_phasediagrams(pids_int)
# [postprocessing(pid) for pid in pids_pd]
# rsync -avz ferpe@atlas-fdr-login-01.sw.ehu.es:/dipc/ferpe/Projects/LMC /Users/fernandopenaranda/Documents/Work/PostdocDonosti/Packages/LMC/src/Paper/ClusterDATA

Figure 9


pdPID = 3238850
drudePID = 3285558
lmcPID = 3286797
ahePID = 3291639
spinPID =3321232

pdpath = "/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Packages/LMC/src/Paper/ClusterDATA/LMC/$(pdPID)_merged_data.jld"
pdpresetpath = "/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Packages/LMC/src/Paper/ClusterDATA/LMC/$(pdPID)_merged_presets.jld"
lmcpath = "/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Packages/LMC/src/Paper/ClusterDATA/LMC/$(lmcPID)_merged_data.jld"
lmcpresetpath = "/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Packages/LMC/src/Paper/ClusterDATA/LMC/$(lmcPID)_merged_presets.jld"
spinpath = "/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Packages/LMC/src/Paper/ClusterDATA/LMC/$(spinPID)_merged_data.jld"
spinpresetpath = "/Users/fernandopenaranda/Documents/Work/PostdocDonosti/Packages/LMC/src/Paper/ClusterDATA/LMC/$(spinPID)_merged_presets.jld"