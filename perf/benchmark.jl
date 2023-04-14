pkgdir = joinpath(@__DIR__, "..")
pkgdir in LOAD_PATH || push!(LOAD_PATH, pkgdir)

import ClimaShallowWater as CSW
using BenchmarkTools
integrator = CSW.setup_integrator(["--output-nsteps", "0"])
trial = BenchmarkTools.@benchmark CSW.step!($integrator)
show(stdout, MIME("text/plain"), trial)
