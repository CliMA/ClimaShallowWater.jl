import ClimaCore.Operators as Operators
import ClimaShallowWater as CSW

using CUDA
using Test
using BenchmarkTools

#####
##### Naive benchmark
#####

function do_work_naive!(y, x)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x
    for i in index:stride:length(y)
        @inbounds y[i] += x[i]
    end
    return
end

function naive_kernel_bm(s)
    # Taken from: https://cuda.juliagpu.org/stable/tutorials/introduction/
    x_d = CUDA.fill(1.0f0, s)
    y_d = CUDA.fill(2.0f0, s)
    N = prod(s)
    N = length(y_d)
    fill!(x_d, 1)
    fill!(y_d, 2)
    kernel = @cuda launch = false do_work_naive!(y_d, x_d)
    config = launch_configuration(kernel.fun)
    threads = min(N, config.threads)
    blocks = cld(N, threads)
    kernel(y_d, x_d; threads, blocks)
    @test all(Array(y_d) .== 3.0f0) # compile and confirm correctness

    # Perform benchmark
    trial = BenchmarkTools.@benchmark $kernel(
        $y_d,
        $x_d;
        threads = $threads,
        blocks = $blocks,
    )
    show(stdout, MIME("text/plain"), trial)
    println()
    return nothing
end

#####
##### ClimaCore benchmark
#####

function do_work_climacore!(integrator, dYdt)
    y = integrator.u
    wdiv = Operators.WeakDivergence()
    grad = Operators.Gradient()
    # Compute Laplacians of h and u
    @. dYdt.h = wdiv(grad(y.h))
    return
end

function climacore_kernel_bm()
    integrator = CSW.setup_integrator(["--output-nsteps", "0"])
    dYdt = copy(integrator.u)
    do_work_climacore!(integrator, dYdt) # compile first, correctness tested in ClimaCore.
    trial = BenchmarkTools.@benchmark do_work_climacore!($integrator, $dYdt)
    show(stdout, MIME("text/plain"), trial)
    println()
    return integrator
end

#####
##### Compare timings
#####

function main()
    integrator = climacore_kernel_bm()
    s = size(parent(integrator.u.h))
    naive_kernel_bm(s)
    return nothing
end
main()

