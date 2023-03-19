
using Makie, CairoMakie, ClimaCoreMakie

function create_movie(
    filename::AbstractString,
    output_dir::AbstractString,
)

    files = sort(filter(endswith(".hdf5"), readdir(output_dir)))

    file = joinpath(output_dir, files[end])
    reader = InputOutput.HDF5Reader(file)
    h = InputOutput.read_field(reader, "Y/h")
    t = InputOutput.HDF5.read_attribute(reader.file, "time")
    Base.close(reader)

    obs = Observable(h)
    scene = plot(obs, axis = (show_axis = false,))

    record(scene, filename, files) do file
        file = joinpath(output_dir, file)
        reader = InputOutput.HDF5Reader(file)
        h = InputOutput.read_field(reader, "Y/h")
        Base.close(reader)
        obs[] = h
    end
end

function create_movie(ARGS::Vector{String}=ARGS)
    s = ArgParseSettings()
    add_arg_group!(s, "Visualization")
    @add_arg_table! s begin
        "--output-dir"
        help = "output directory"
        default = "output"
        "--filename"
        help = "output filename"
        default = "output/out.mp4"
    end
    parsed_args = parse_args(ARGS, s)
    create_movie(
        parsed_args["filename"],
        parsed_args["output-dir"], 
    )
end