
using ArgParse

function setup_integrator(
    space::Spaces.SpectralElementSpace2D, 
    test::AbstractSphereTestCase; 
    time_step=60*6, 
    time_end=60*60*2,
    output_nsteps=5,
    output_dir="output",
    )

    f = coriolis_field(space, test)
    h_s = surface_height_field(space, test)
    Y = initial_condition(space, test)
    p = auxiliary_state(Y, test)
    dss!(Y, p) # ensure continuous

    FT = Spaces.undertype(space)

    # Solve the ODE
    if save_nsteps > 0
        save_callback = ClimaTimeSteppers.Callbacks.EveryXSimulationSteps(
            save_nsteps;
            atinit = true,
        ) do integrator
            NVTX.@range "saving output" begin
                output_dir = "output"
                Y = integrator.u
                t = integrator.t
                day = floor(Int, t / (60 * 60 * 24))
                sec = floor(Int, t % (60 * 60 * 24))
                @info "Saving state" day sec
                mkpath(output_dir)
                output_file = joinpath(output_dir, "day$day.$sec.hdf5")
                hdfwriter = InputOutput.HDF5Writer(output_file, space.topology.context)
                InputOutput.HDF5.write_attribute(hdfwriter.file, "time", t)
                InputOutput.write!(hdfwriter, "Y" => Y)
                Base.close(hdfwriter)
            end
            return nothing
        end
    else
        save_callback = nothing
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
        callback = CallbackSet(save_callback),
    )

    return integrator
end


function setup_integrator(ARGS::Vector{String}=ARGS)
    CUDA.allowscalar(false)

    s = ArgParseSettings()

    @add_arg_table! s begin
        "--float-type"
            help = "Floating point type to use"
            eval_arg = true
            default = Float64
        "--panel-size"
            help = "Number of elements across each panel"
            arg_type = Integer
            default = 8
        "--poly-nodes"
            help = "Number of nodes in each dimension to use in the polynomial approximation. Polynomial degree = poly-nodes - 1."
            arg_type = Integer
            default = 4
        "--time-step"
            help = "Time step (seconds)"
            arg_type = Float64
            default = Float64(60*6)
        "--time-end"
            help = "End time (seconds)"
            arg_type = Float64
            default = Float64(60*60*24*2)
        "--save-nsteps"
            help = "Number of time steps between saving output"
            arg_type = Integer
            default = 5
        "testcase"
            help = "Test case to run"
            eval_arg = true
            default = SteadyStateTest()
    end
    args = parse_args(ARGS, s)

    testcase = args["testcase"]

    space = create_space(
        testcase;
        float_type = args["float-type"],
        panel_size = args["panel-size"],
        poly_nodes = args["poly-nodes"],
    )

    setup_integrator(
        space, 
        testcase; 
        time_step = args["time-step"],
        time_end = args["time-end"],
        save_nsteps = args["save-nsteps"],
    )
end


function run(ARGS::Vector{String}=ARGS)
    integrator = setup_integrator(ARGS)
    solve!(integrator)
end


