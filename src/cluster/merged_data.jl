#using LMC, JLD2 # needs to load LMC structures to approately read the preset files
""" once the slurm calculation if finished we run these julia functions on the cluster to create merged data files"""
function data_merge(folderPID)
    merged = Dict{String,Any}()
    pathtoPID = find_folder(folderPID) # run in LMC folder project creates the absolute path to PID
    pathtofiles = glob(joinpath("*", "data.jld"), pathtoPID) 
    sorted_paths = sort( pathtofiles, by = p -> parse(Int, basename(dirname(p)))) # sort numerically by folder order which is Ez order
    pathtopresets = glob(joinpath("*", "presets.jld"), pathtoPID) 
    for f in sorted_paths
        merged[f] = load(f)
    end
    presets = load(pathtopresets[1])
    @save pathtoPID * "/merged_pd_data.jld" merged # data
    @save pathtoPID * "/merged_presets.jld" presets
end

""" PID folder finder: e.g. PID = 0000001, 
    `find_folder(PID)`
 find the absolute path to PID folder regardless the data
 store phase diagrams, drude, qah, LMCS responses... """
find_folder(PID) = find_folder(string(PID))
function find_folder(target::AbstractString)
    root = pwd()
    for (dirpath, dirnames, _) in walkdir(root)
        if target in dirnames
            return joinpath(dirpath, target)
        end
    end
    return nothing
end