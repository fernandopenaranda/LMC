include("bistritzermodel.jl")
include("non_linear_responses.jl")
include("data_processing.jl")
include("plots_bistritzer.jl")
include("tests.jl")
include("k_resolved_calculations.jl")

function plot2valleybands(p, mass_term, ylimits = [-0.06,0.06]; resolution = (400,300), ph = false)
    es = bands_bistritzer(ParamsBM(p , ν = 1, pointsk = 50), eigvecs = false, mass_term = mass_term, ph = ph)[1];
    fig = plotbands(es, color = :gray, ylimits = ylimits, resolution = resolution)
    es = bands_bistritzer(ParamsBM(p , ν = -1, pointsk = 50), eigvecs = false, mass_term = mass_term, ph = ph)[1];
    # plotbands!(fig, es, color = :orange, ylimits = ylimits)
    return fig
end