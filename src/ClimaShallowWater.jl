module ClimaShallowWater

using ClimaComms, ClimaCommsMPI
using DocStringExtensions
using LinearAlgebra
using ClimaTimeSteppers, DiffEqBase
using DiffEqCallbacks
import CUDA
using NVTX
using ArgParse


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
