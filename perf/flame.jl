import Profile
import ClimaShallowWater as CSW

function do_work!(integrator)
    for _ in 1:100
        CSW.step!(integrator)
    end
    return nothing
end

integrator = CSW.setup_integrator(["--output-nsteps", "0"])
do_work!(integrator) # compile first
Profile.clear_malloc_data()
prof = Profile.@profile begin
    do_work!(integrator)
end

import ProfileCanvas

if haskey(ENV, "BUILDKITE_COMMIT") || haskey(ENV, "BUILDKITE_BRANCH")
    output_dir = joinpath(@__DIR__, "output")
    mkpath(output_dir)
    ProfileCanvas.html_file(joinpath(output_dir, "flame.html"))
else
    ProfileCanvas.view(Profile.fetch())
end
