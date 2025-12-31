#_________________________________________________________________________________________
# LMC
#_________________________________________________________________________________________
"""
function observable functions are defined for each flavour so you have to sum
the four of them to get the experimental value
"""
function lmc_map(N, p::Params_rhombohedral, Ezlist, νlist, μαs; kws...)
    presets(; ξ, Ez, μ) = xxx_lmc_presets(N, 
        Params_rhombohedral(p, ξ = ξ, μ = μ, Delta_Ez = Ez); kws...)
    return map_eval(linear_magneto_conductivity_orbital, presets, μαs, Ezlist, νlist)
end

#_______________________________________________________________________________________
# Drude
#_________________________________________________________________________________________
function drude_map(N, p::Params_rhombohedral, Ezlist, νlist, μαs; kws...)
    presets(; ξ, Ez, μ) = 
        xx_drude_presets(N, Params_rhombohedral(p, ξ = ξ, μ = μ, Delta_Ez = Ez); kws...)
        return map_eval(drude_conductivity, presets, μαs, Ezlist, νlist)
end


#_________________________________________________________________________________________
# AHE. xy
#_________________________________________________________________________________________
function ahe_map(N, p::Params_rhombohedral, Ezlist, νlist, μαs; kws...)
    presets(; ξ, Ez, μ) = 
        qah_presets(N, :x, :y, Params_rhombohedral(p, ξ = ξ, μ = μ, Delta_Ez = Ez); kws...)
    return map_eval(σij_anomalous_hall, presets, μαs, Ezlist, νlist)
end

#_________________________________________________________________________________________
# Evaluator
#_________________________________________________________________________________________
function map_eval(obs, presets::Union{Planar_σijk_presets_orbital, AH_presets, σij_presets}, μαs, Ezlist, νlist)
    μs = reshape_densities(μαs)
    dimx = size(μs[1],1)
    dimy = size(μs[1],2)
    ξs = [1, -1, 1, -1] # valley array
    mat = zeros(Float64, dimx, dimy)
    for (i,Ez) in enumerate(Ezlist)
        for (j, ν) in enumerate(νlist)
            mat[i,j] = sum([obs(presets(ξ = ξs[k], Ez = Ez, μ = μs[k][i,j])) for k in 1:4]) # sum over the 4 flavors
            # mat[i,j] = obs(presets(ξ = ξs[1], Ez = Ez, μ = μs[1][i,j])) # single flavour
        end
    end
    return mat
end
#separate method for spin contribution to the lmc
function map_eval(obs, presets::Planar_σijk_presets_spin, μαs, Ezlist, νlist)
    μs = reshape_densities(μαs)
    dimx = size(μs[1],1)
    dimy = size(μs[1],2)
    ξs = [1, -1, 1, -1] # valley array
    σs = [1, 1, -1, -1] # spin array
    mat = zeros(Float64, dimx, dimy)
    for (i,Ez) in enumerate(Ezlist)
        for (j, ν) in enumerate(νlist)
            mat[i,j] = sum([σs[k] * obs(presets(ξ = ξs[k], Ez = Ez, μ = μs[k][i,j])) for k in 1:4]) # sum over the 4 flavors note the spin prefactor due to Hunds coupling
        end
    end
    return mat
end

function Ez_map_eval(obs, presets, μs, Ez, νlist)
    ξs = [1, -1, 1, -1]
    obs_list = zeros(length(νlist))
    if obs == linear_magneto_conductivity_spin
        σs = [1, 1, -1, -1] # spin array
        for (j, ν) in enumerate(νlist)
            obs_list[j] = sum([σs[k] * obs(presets(ξ = ξs[k], Ez = Ez, μ = μs[j][k]))
                 for k in 1:4]) # sum over the 4 flavors
         end
    else
        for (j, ν) in enumerate(νlist)
            obs_list[j] = sum([obs(presets(ξ = ξs[k], Ez = Ez, μ = μs[j][k])) 
                 for k in 1:4]) # sum over the 4 flavors
         end
    end
    return obs_list
end
