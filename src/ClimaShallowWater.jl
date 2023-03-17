module ClimaShallowWater

using ClimaComms
using DocStringExtensions
using LinearAlgebra
using ClimaTimeSteppers, DiffEqBase
using DiffEqCallbacks
import CUDA
using NVTX

import ClimaCore:
    Device,
    Domains,
    Fields,
    Geometry,
    Meshes,
    Operators,
    Spaces,
    Topologies,
    DataLayouts,
    InputOutput


include("common.jl")
include("testcases.jl")
include("driver.jl")

end # module ClimaShallowWater
