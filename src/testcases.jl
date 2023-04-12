import CLIMAParameters as CP

function construct_paramset(
    testcase::Type{T},
    toml_dict,
) where {T <: AbstractSphereTestCase}
    param_names = string.(fieldnames(testcase))
    parameters = CP.get_parameter_values!(toml_dict, param_names)
    FT = CP.float_type(toml_dict)
    spherical_parameters = SphericalParameters(toml_dict)
    T{FT}(; parameters..., common_parameters = spherical_parameters)
end

"""
    SphericalParameters

Physical parameters needed for all spherical simulations

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SphericalParameters{FT}
    "Radius of earth"
    planet_radius::FT
    "Rotation rate of earth"
    rotation_rate::FT
    "Gravitational constant"
    grav::FT
    "Hyperdiffusion coefficient"
    hyperdiff_coefficient::FT
    "angle between the north pole and the center of the top cube panel"
    angle_α::FT
end
function SphericalParameters(toml_dict)
    param_names = string.(fieldnames(SphericalParameters))
    parameters = CP.get_parameter_values!(toml_dict, param_names)
    FT = CP.float_type(toml_dict)
    SphericalParameters{FT}(; parameters...)
end

"""
    SteadyStateTest()

Test case 2 of [Williamson1992](@cite).

 This test case gives the steady-state solution to the non-linear shallow water
equations. It consists of solid body rotation or zonal flow with the
corresponding geostrophic height field. This can be run with an angle α that
represents the angle between the north pole and the center of the top cube
panel.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SteadyStateTest{FT} <: AbstractSphereTestCase
    "Physical parameters"
    common_parameters::SphericalParameters{FT}
    "advection velocity"
    advection_velocity::FT = 2 * pi * common_parameters.planet_radius / (12 * 86400)
    "constant for computing peak of analytic height field"
    peak_analytic_height_field_parameter::FT
    "peak of analytic height field"
    peak_analytic_height_field::FT = peak_analytic_height_field_parameter / common_parameters.grav
end

function initial_height(space, test::SteadyStateTest)
    u0 = advection_velocity(test)
    h0 = peak_analytic_height_field(test)
    R = planet_radius(test)
    Ω = rotation_rate(test)
    g = grav(test)
    α = angle_α(test)

    coordinates = Fields.coordinate_field(space)
    ϕ = coordinates.lat
    λ = coordinates.long
    h = @. h0 -
       (R * Ω * u0 + u0^2 / 2) / g *
       (-cosd(λ) * cosd(ϕ) * sind(α) + sind(ϕ) * cosd(α))^2
    return h
end
function initial_velocity(space, test::SteadyStateTest)
    u0 = advection_velocity(test)
    α = angle_α(test)

    coordinates = Fields.coordinate_field(space)
    ϕ = coordinates.lat
    λ = coordinates.long

    uλ = @. u0 * (cosd(α) * cosd(ϕ) + sind(α) * cosd(λ) * sind(ϕ))
    uϕ = @. -u0 * sind(α) * sind(λ)

    u = @. Geometry.Covariant12Vector(Geometry.UVVector(uλ, uϕ))
    return u
end



"""
    SteadyStateCompactTest()

Test case 3 of [Williamson1992](@cite).

This test case gives the steady-state solution to the
non-linear shallow water equations with nonlinear zonal geostrophic flow
with compact support.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SteadyStateCompactTest{FT} <: AbstractSphereTestCase
    "Physical parameters"
    common_parameters::SphericalParameters{FT}
    "latitude lower bound for coordinate transformation parameter"
    coord_transform_lower_bound::FT
    "latitude upper bound for coordinate transformation parameter"
    coord_transform_upper_bound::FT
    "velocity perturbation parameter"
    velocity_perturbation::FT
    "advection velocity"
    advection_velocity::FT = 2 * pi * common_parameters.planet_radius / (12 * 86400)
    "constant for computing peak of analytic height field"
    peak_analytic_height_field_parameter::FT
    "peak of analytic height field"
    peak_analytic_height_field::FT = peak_analytic_height_field_parameter / common_parameters.grav
end

function initial_condition(
    space::Spaces.SpectralElementSpace2D,
    test::SteadyStateCompactTest,
)
    u0 = advection_velocity(test)
    h0 = peak_analytic_height_field(test)
    ϕᵦ = coord_transform_lower_bound(test)
    ϕₑ = coord_transform_upper_bound(test)
    xₑ = velocity_perturbation(test)
    R = planet_radius(test)
    Ω = rotation_rate(test)
    g = grav(test)
    α = angle_α(test)

    Y = map(Fields.local_geometry_field(space)) do local_geometry
        coord = local_geometry.coordinates

        ϕ = coord.lat
        λ = coord.long

        if α == 0.0
            ϕprime = ϕ
            λprime = λ
        else
            ϕprime = asind(sind(ϕ) * cosd(α) - cosd(ϕ) * cosd(λ) * sind(α))
            λprime = asind(sind(λ) * cosd(ϕ) / cosd(ϕprime)) # for alpha45, this experiences numerical precision issues. The test case is designed for either alpha0 or alpha60

            # Temporary angle to ensure λprime is in the right quadrant
            λcond = cosd(α) * cosd(λ) * cosd(ϕ) + sind(α) * sind(ϕ)

            # If λprime is not in the right quadrant, adjust
            if λcond < 0.0
                λprime = -λprime - 180.0 # shifted by 180 compared to the paper, because our λ ∈ [-180, 180]
            end
            if λprime < -180.0
                λprime += 360.0
            end
        end

        # Set auxiliary function needed for initial state of velocity field
        b(x) = x ≤ 0.0 ? 0.0 : exp(-x^(-1))

        x(ϕprime) = xₑ * (ϕprime - ϕᵦ) / (ϕₑ - ϕᵦ)
        uλprime(ϕprime) =
            u0 * b(x(ϕprime)) * b(xₑ - x(ϕprime)) * exp(4.0 / xₑ)
        uϕprime = 0.0

        # Set integral needed for height initial state
        h_int(γ) =
            abs(γ) < 90.0 ?
            (2 * Ω * sind(γ) + uλprime(γ) * tand(γ) / R) * uλprime(γ) : 0.0

        # Set initial state for height field
        h =
            h0 - (R / g) * (pi / 180.0) * QuadGK.quadgk(h_int, -90.0, ϕprime)[1]

        # Set initial state for velocity field
        uϕ = -(uλprime(ϕprime) * sind(α) * sind(λprime)) / cosd(ϕ)
        if abs(cosd(λ)) < 1e-13
            if abs(α) > 1e-13
                if cosd(λ) > 0.0
                    uλ = -uϕ * cosd(ϕ) / tand(α)
                else
                    uλ = uϕ * cosd(ϕ) / tand(α)
                end
            else
                uλ = uλprime(ϕprime)
            end
        else
            uλ =
                (uϕ * sind(ϕ) * sind(λ) + uλprime(ϕprime) * cosd(λprime)) /
                cosd(λ)
        end

        u = Geometry.transform(
            Geometry.Covariant12Axis(),
            Geometry.UVVector(uλ, uϕ),
            local_geometry,
        )

        return (h = h, u = u)
    end
    return Y
end

"""
    MountainTest()

Test case 5 of [Williamson1992](@cite).

It represents a zonal flow over an isolated mountain,
where the governing equations describe a global steady-state nonlinear
zonal geostrophic flow, with a corresponding geostrophic height field over
a non-uniform reference surface h_s.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct MountainTest{FT} <: AbstractSphereTestCase
    "Physical parameters"
    common_parameters::SphericalParameters{FT}
    "advection velocity"
    advection_velocity::FT
    "peak of analytic height field"
    peak_analytic_height_field::FT
    "radius of conical mountain"
    velocity_amplitude::FT
    "center of mountain long coord, shifted by 180 compared to the paper, 
    because our ``λ ∈ [-180, 180]`` (in the paper it was 270, with ``λ ∈ [0, 360]``)"
    mountain_center_longitude::FT
    "latitude coordinate for center of mountain"
    mountain_center_latitude::FT
    "mountain peak height"
    mountain_peak_height::FT
end

function surface_height_field(
    space::Spaces.SpectralElementSpace2D,
    test::MountainTest,
)
    a = velocity_amplitude(test)
    λc = mountain_center_longitude(test)
    ϕc = mountain_center_latitude(test)
    h_s0 = mountain_peak_height(test)
    map(Fields.coordinate_field(space)) do coord
        ϕ = coord.lat
        λ = coord.long
        r = sqrt(min(a^2, (λ - λc)^2 + (ϕ - ϕc)^2)) # positive branch
        h_s0 * (1 - r / a)
    end
end
function initial_condition(
    space::Spaces.SpectralElementSpace2D,
    test::MountainTest{FT},
)   where {FT}
    # steady-state and mountain test cases share the same initial condition
    peak_analytic_height_field_parameter = peak_analytic_height_field(test) * grav(test)
    steady_state_test = SteadyStateTest{FT}(
        common_parameters(test),
        advection_velocity(test),
        peak_analytic_height_field_parameter,
        peak_analytic_height_field(test)
    )
    initial_condition(space, steady_state_test)
end


"""
    RossbyHaurwitzTest <: AbstractSphereTestCase

Test case 6 of [Williamson1992](@cite).

It represents the solution of the nonlinear barotropic
vorticity equation on the sphere

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct RossbyHaurwitzTest{FT} <: AbstractSphereTestCase
    "Physical parameters"
    common_parameters::SphericalParameters{FT}
    "velocity amplitude parameter"
    velocity_amplitude::FT
    "peak of analytic height field"
    peak_analytic_height_field::FT
    "vorticity amplitude parameter (1/sec)"
    vorticity_amplitude_ω::FT
    "vorticity amplitude parameter (1/sec)"
    vorticity_amplitude_K::FT
end

function initial_height(
    space::Spaces.SpectralElementSpace2D,
    test::RossbyHaurwitzTest,
)
    a = velocity_amplitude(test)
    h0 = peak_analytic_height_field(test)
    ω = vorticity_amplitude_ω(test)
    K = vorticity_amplitude_K(test)
    R = planet_radius(test)
    Ω = rotation_rate(test)
    g = grav(test)

    map(Fields.local_geometry_field(space)) do local_geometry
        coord = local_geometry.coordinates
        ϕ = coord.lat
        λ = coord.long

        A =
            ω / 2 * (2 * Ω + ω) * cosd(ϕ)^2 +
            1 / 4 *
            K^2 *
            (
                (a + 1) * cosd(ϕ)^(2 * a + 2) +
                (2 * a^2 - a - 2) * cosd(ϕ)^(2 * a) -
                2 * a^2 * cosd(ϕ)^(2 * a - 2)
            )
        B =
            2 * (Ω + ω) * K / (a + 1) / (a + 2) *
            cosd(ϕ)^a *
            ((a^2 + 2 * a + 2) - (a + 1)^2 * cosd(ϕ)^2)
        C = 1 / 4 * K^2 * cosd(ϕ)^(2 * a) * ((a + 1) * cosd(ϕ)^2 - (a + 2))

        h =
            h0 +
            (R^2 * A + R^2 * B * cosd(a * λ) + R^2 * C * cosd(2 * a * λ)) / g
    end
end
function initial_velocity(
    space::Spaces.SpectralElementSpace2D,
    test::RossbyHaurwitzTest,
)
    a = velocity_amplitude(test)
    h0 = peak_analytic_height_field(test)
    ω = vorticity_amplitude_ω(test)
    K = vorticity_amplitude_K(test)
    R = planet_radius(test)
    Ω = rotation_rate(test)
    g = grav(test)

    map(Fields.local_geometry_field(space)) do local_geometry
        coord = local_geometry.coordinates
        ϕ = coord.lat
        λ = coord.long


        uλ =
            R * ω * cosd(ϕ) +
            R * K * cosd(ϕ)^(a - 1) * (a * sind(ϕ)^2 - cosd(ϕ)^2) * cosd(a * λ)
        uϕ = -R * K * a * cosd(ϕ)^(a - 1) * sind(ϕ) * sind(a * λ)


        Geometry.transform(
            Geometry.Covariant12Axis(),
            Geometry.UVVector(uλ, uϕ),
            local_geometry,
        )
    end
end



"""
    BarotropicInstabilityTest{FT} <: AbstractSphereTestCase

Test case from [Galewsky2004](@cite) (and also §7.6 of [Ullrich2010](@cite)).

This test case consists of a zonal jet with compact
support at a latitude of 45°. A small height disturbance is then added,
which causes the jet to become unstable and collapse into a highly vortical
structure.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct BarotropicInstabilityTest{FT} <: AbstractSphereTestCase
    "Physical parameters"
    common_parameters::SphericalParameters{FT}
    "maximum zonal velocity"
    max_zonal_velocity::FT
    "mountain shape parameters"
    mountain_shape_α::FT
    "mountain shape parameters"
    mountain_shape_β::FT
    "peak of balanced height field from Tempest 
    https://github.com/paullric/tempestmodel/blob/master/test/shallowwater_sphere/BarotropicInstabilityTest.cpp#L86"
    peak_analytic_height_field::FT
    "local perturbation peak height"
    local_perturb_peak_height::FT
    "southern jet boundary"
    southern_jet_boundary::FT
    "northern jet boundary"
    northern_jet_boundary::FT
    "height perturbation peak location"
    height_perturbation_peak_location::FT
end

"zonal velocity decay parameter"
zonal_velocity_decay(test::BarotropicInstabilityTest) = exp(-4.0 / (deg2rad(northern_jet_boundary(test)) - deg2rad(southern_jet_boundary(test)))^2)

function initial_height(
    space::Spaces.SpectralElementSpace2D,
    test::BarotropicInstabilityTest,
)
    # we need to instantiate the initial height field on the CPU so that we can use quadgk
    initial_height(space, test, Device.device(space))
end

function initial_height(
    space::Spaces.SpectralElementSpace2D,
    test::BarotropicInstabilityTest,
    ::ClimaComms.CUDA,
)
    mesh = space.topology.mesh
    context = space.topology.context
    if context isa ClimaComms.SingletonCommsContext
        cpu_context = ClimaComms.SingletonCommsContext(ClimaComms.CPU())
    elseif context isa ClimaCommsMPI.MPICommsContext
        cpu_context = ClimaComms.MPICommsContext(ClimaComms.CPU())
    end
    cpu_topology = Topologies.Topology2D(cpu_context, mesh)
    cpu_space = Spaces.SpectralElementSpace2D(
        cpu_topology,
        Spaces.quadrature_style(space),
    )
    cpu_h = initial_height(cpu_space, test)
    h = zeros(space)
    copyto!(parent(h), parent(cpu_h))
    return h
end



function initial_height(
    space::Spaces.SpectralElementSpace2D,
    test::BarotropicInstabilityTest,
    ::ClimaComms.CPU,
)
    u_max = max_zonal_velocity(test)
    αₚ = mountain_shape_α(test)
    βₚ = mountain_shape_β(test)
    h0 = peak_analytic_height_field(test)
    h_hat = local_perturb_peak_height(test)
    ϕ₀ = southern_jet_boundary(test)
    ϕ₁ = northern_jet_boundary(test)
    ϕ₂ = height_perturbation_peak_location(test)
    eₙ = zonal_velocity_decay(test)
    R = planet_radius(test)
    Ω = rotation_rate(test)
    g = grav(test)
    α = angle_α(test)

    map(Fields.local_geometry_field(space)) do local_geometry
        coord = local_geometry.coordinates

        ϕ = coord.lat
        λ = coord.long

        if α == 0.0
            ϕprime = ϕ
        else
            ϕprime = asind(sind(ϕ) * cosd(α) - cosd(ϕ) * cosd(λ) * sind(α))
        end

        # Set initial state of velocity field
        uλprime(ϕprime) =
            (u_max / eₙ) *
            exp(1.0 / (deg2rad(ϕprime - ϕ₀) * deg2rad(ϕprime - ϕ₁))) *
            (ϕ₀ < ϕprime < ϕ₁)
        uϕprime = 0.0

        # Set integral needed for height initial state
        h_int(γ) =
            abs(γ) < 90.0 ?
            (2 * Ω * sind(γ) + uλprime(γ) * tand(γ) / R) * uλprime(γ) : 0.0

        # Set initial state for height field
        h = h0 - (R / g) * (pi / 180.0) * QuadGK.quadgk(h_int, -90.0, ϕprime)[1]

        if λ > 0.0
            λ -= 360.0
        end
        if λ < -360.0 || λ > 0.0
            @info "Invalid longitude value"
        end

        # Add height perturbation
        h += h_hat * cosd(ϕ) * exp(-(λ^2 / αₚ^2) - ((ϕ₂ - ϕ)^2 / βₚ^2))

        return h
    end
end

function initial_velocity(
    space::Spaces.SpectralElementSpace2D,
    test::BarotropicInstabilityTest,
)
    u_max = max_zonal_velocity(test)
    αₚ = mountain_shape_α(test)
    βₚ = mountain_shape_β(test)
    h0 = peak_analytic_height_field(test)
    h_hat = local_perturb_peak_height(test)
    ϕ₀ = southern_jet_boundary(test)
    ϕ₁ = northern_jet_boundary(test)
    ϕ₂ = height_perturbation_peak_location(test)
    eₙ = zonal_velocity_decay(test)
    R = planet_radius(test)
    Ω = rotation_rate(test)
    g = grav(test)
    α = angle_α(test)

    map(Fields.local_geometry_field(space)) do local_geometry
        coord = local_geometry.coordinates

        ϕ = coord.lat
        λ = coord.long

        if α == 0.0
            ϕprime = ϕ
        else
            ϕprime = asind(sind(ϕ) * cosd(α) - cosd(ϕ) * cosd(λ) * sind(α))
        end

        # Set initial state of velocity field
        uλprime(ϕprime) =
            (u_max / eₙ) *
            exp(1.0 / (deg2rad(ϕprime - ϕ₀) * deg2rad(ϕprime - ϕ₁))) *
            (ϕ₀ < ϕprime < ϕ₁)
        uϕprime = 0.0

        uλ = uλprime(ϕprime)
        uϕ = uϕprime

        u = Geometry.transform(
            Geometry.Covariant12Axis(),
            Geometry.UVVector(uλ, uϕ),
            local_geometry,
        )
    end
end

# Generate functions for accessing parameters
for fn in fieldnames(SphericalParameters)
    @eval $(fn)(ps::SphericalParameters) = ps.$(fn)
end

common_params(ps::AbstractSphereTestCase) = ps.common_parameters

for paramset in (
    SteadyStateTest,
    SteadyStateCompactTest,
    MountainTest,
    RossbyHaurwitzTest,
    BarotropicInstabilityTest,
)
    for var in fieldnames(paramset)
        @eval $(var)(ps::$(paramset)) = ps.$(var)
    end
    for var in fieldnames(SphericalParameters)
        @eval $var(ps::$(paramset)) = $var(common_params(ps))
    end
end
