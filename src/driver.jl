function setup_integrator(
    space::Spaces.SpectralElementSpace2D, 
    test::AbstractSphereTestCase; 
    time_step=60*6, 
    time_end=60*60*2,
    output_nsteps=5,
    output_dir="output",
    )

    Y = initial_condition(space, test)
    p = auxiliary_state(Y, test)
    dss!(Y, p) # ensure continuous

    FT = Spaces.undertype(space)

    context = space.topology.context

    # Solve the ODE
    if output_nsteps > 0
        mkpath(output_dir)
        output_n = 0
        output_callback = ClimaTimeSteppers.Callbacks.EveryXSimulationSteps(
            output_nsteps;
            atinit = true,
        ) do integrator
            global n
            NVTX.@range "saving output" begin
                Y = integrator.u
                t = integrator.t
                output_file = joinpath(output_dir, "state_$(string(output_n, pad=6)).hdf5")
                if ClimaComms.iamroot(context)
                    @info "Saving state" n=output_n output_file t
                end            
                hdfwriter = InputOutput.HDF5Writer(output_file, space.topology.context)
                InputOutput.HDF5.write_attribute(hdfwriter.file, "time", t)
                InputOutput.write!(hdfwriter, "Y" => Y)
                Base.close(hdfwriter)
                output_n += 1
            end
            return nothing
        end
    else
        output_callback = nothing
    end

    prob =
        ODEProblem(ClimaODEFunction(; T_exp! = tendency!), Y, (FT(0), FT(time_end)), p)
    integrator = DiffEqBase.init(
        prob,
        ExplicitAlgorithm(SSP33ShuOsher()),
        dt = FT(time_step),
        saveat = [],
        progress = true,
        adaptive = false,
        progress_message = (dt, u, p, t) -> t,
        internalnorm = (u, t) -> norm(parent(Y)),
        callback = CallbackSet(output_callback),
    )

    return integrator
end


function ismpi()
    # detect common environment variables used by MPI launchers
    #   PMI_RANK appears to be used by MPICH and srun
    #   OMPI_COMM_WORLD_RANK appears to be used by OpenMPI
    return haskey(ENV, "PMI_RANK") || haskey(ENV, "OMPI_COMM_WORLD_RANK")
end
function setup_integrator(ARGS::Vector{String}=ARGS)
    CUDA.allowscalar(false)

    s = ArgParseSettings(prog="shallowwater")

    @add_arg_table! s begin
        "--device"
            help = "Computation device (CPU, CUDA)"
            arg_type = String
            default = CUDA.functional() ? "CUDA" : "CPU"
        "--comms"
            help = "Communication type (Singleton, MPI)"
            arg_type = String
            default = ismpi() ? "MPI" : "Singleton"
        "--float-type"
            help = "Floating point type (Float32, Float64)"
            eval_arg = true
            default = Float64
        "--panel-size"
            help = "Number of elements across each panel"
            arg_type = Int
            default = 8
        "--poly-nodes"
            help = "Number of nodes in each dimension to use in the polynomial approximation. Polynomial degree = poly-nodes - 1."
            arg_type = Int
            default = 4
        "--time-step"
            help = "Time step (seconds)"
            arg_type = Float64
            default = Float64(60*6)
        "--time-end"
            help = "End time (seconds)"
            arg_type = Float64
            default = Float64(60*60*24*2)
        "--output-nsteps"
            help = "Number of time steps between saving output"
            arg_type = Integer
            default = 5
        "--output-dir"
            help = "Directory to save output to"
            arg_type = String
            default = "output"
        "testcase"
            help = "Test case to run"
            default = "steadystate"
    end
    args = parse_args(ARGS, s)

    device = args["device"] == "CUDA" ? ClimaComms.CUDA() : 
             args["device"] == "CPU" ? ClimaComms.CPU() :
             error("Unknown device: $(args["device"])")

    context = args["comms"] == "MPI" ? ClimaCommsMPI.MPICommsContext(device) :
              args["comms"] == "Singleton" ? ClimaComms.SingletonCommsContext(device) :
              error("Unknown comms: $(args["comms"])")

    ClimaComms.init(context)

    testcase = args["testcase"] == "steadystate" ? SteadyStateTest() :
    args["testcase"] == "steadystatecompact" ? SteadyStateCompactTest() :
        args["testcase"] == "mountain" ? MountainTest() :
        args["testcase"] == "rossbyhaurwitz" ? RossbyHaurwitzTest() :
        args["testcase"] == "baroclinic" ? BaroclinicWaveTest() :
    error("Unknown testcase: $(args["testcase"])")
    float_type = args["float-type"]
    panel_size = args["panel-size"]
    poly_nodes = args["poly-nodes"]
    time_step = args["time-step"]
    time_end = args["time-end"]


    space = create_space(
        context,
        testcase;
        float_type,
        panel_size,
        poly_nodes,
    )

    if ClimaComms.iamroot(context)
        nprocs = ClimaComms.nprocs(context)
        @info "Setting up experiment" device context testcase float_type panel_size poly_nodes time_step time_end approx_resolution=approx_resolution(space) D₄ = hyperdiffusion_coefficient(space, testcase)
    end
    setup_integrator(
        space,
        testcase; 
        time_step,
        time_end,
        output_nsteps = args["output-nsteps"],
        output_dir = args["output-dir"],
    )
end


function run(ARGS::Vector{String}=ARGS)
    integrator = setup_integrator(ARGS)
    solve!(integrator)
end


