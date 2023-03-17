"""
    SphericalParameters

Physical parameters needed for all spherical simulations

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SphericalParameters
    "Radius of earth"
    R::Float64 = 6.37122e6
    "Rotation rate of earth"
    Ω::Float64 = 7.292e-5
    "Gravitational constant"
    g::Float64 = 9.80616
    "Hyperdiffusion coefficient"
    ν::Float64 = 0.0015
    "angle between the north pole and the center of the top cube panel"
    α::Float64 = 0.0
end


"""
    SteadyStateTest()

The first one, called "steady_state", reproduces Test Case 2 in Williamson et al,
"A standard test set for numerical approximations to the shallow water
equations in spherical geometry", 1992. This test case gives the steady-state
solution to the non-linear shallow water equations. It consists of solid
body rotation or zonal flow with the corresponding geostrophic height field.
This can be run with an angle α that represents the angle between the north
pole and the center of the top cube panel.

https://doi.org/10.1016/S0021-9991(05)80016-6

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SteadyStateTest <: AbstractSphereTestCase
    "Physical parameters"
    params::SphericalParameters = SphericalParameters()
    "advection velocity"
    u0::Float64 = 2 * pi * params.R / (12*60*60*24) # 1 revolution/12 days
    "peak of analytic height field"
    h0::Float64 = 2.94e4 / params.g
end

function initial_height(
    space,
    test::SteadyStateTest,
)
    FT = Spaces.undertype(space)
    u0 = FT(test.u0)
    h0 = FT(test.h0)
    R = FT(test.params.R)
    Ω = FT(test.params.Ω)
    g = FT(test.params.g)
    α = FT(test.params.α)

    coordinates = Fields.coordinate_field(space)
    ϕ = coordinates.lat
    λ = coordinates.long
    h = @. h0 - (R * Ω * u0 + u0^2 / 2) / g *
        (-cosd(λ) * cosd(ϕ) * sind(α) + sind(ϕ) * cosd(α))^2
    return h
end
function initial_velocity(
    space,
    test::SteadyStateTest,
)
    FT = Spaces.undertype(space)
    u0 = FT(test.u0)
    h0 = FT(test.h0)
    R = FT(test.params.R)
    Ω = FT(test.params.Ω)
    g = FT(test.params.g)
    α = FT(test.params.α)

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

A second one, called "steady_state_compact", reproduces Test Case 3 in the same
reference paper. This test case gives the steady-state solution to the
non-linear shallow water equations with nonlinear zonal geostrophic flow
with compact support.

https://doi.org/10.1016/S0021-9991(05)80016-6

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SteadyStateCompactTest{FT, P} <: AbstractSphereTestCase
    "Physical parameters"
    params::SphericalParameters = SphericalParameters()
    "advection velocity"
    u0::Float64 = 2 * pi * params.R / (12 * 86400)
    "peak of analytic height field"
    h0::Float64 = 2.94e4 / params.g
    "latitude lower bound for coordinate transformation parameter"
    ϕᵦ::Float64 = -30.0
    "latitude upper bound for coordinate transformation parameter"
    ϕₑ::Float64 = 90.0
    "velocity perturbation parameter"
    xₑ::Float64 = 0.3
end

function initial_condition(space, test::SteadyStateCompactTest)
    FT = Spaces.undertype(space)
    u0 = FT(test.u0)
    h0 = FT(test.h0)
    ϕᵦ = FT(test.ϕᵦ)
    ϕₑ = FT(test.ϕₑ)
    xₑ = FT(test.xₑ)
    R = FT(test.params.R)
    Ω = FT(test.params.Ω)
    g = FT(test.params.g)
    α = FT(test.params.α)

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

A third one, called "mountain", reproduces Test Case 5 in the same
reference paper. It represents a zonal flow over an isolated mountain,
where the governing equations describe a global steady-state nonlinear
zonal geostrophic flow, with a corresponding geostrophic height field over
a non-uniform reference surface h_s.

https://doi.org/10.1016/S0021-9991(05)80016-6

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct MountainTest <: AbstractSphereTestCase
    "Physical parameters"
    params::SphericalParameters = SphericalParameters()
    "advection velocity"
    u0::Float64 = 20.0
    "peak of analytic height field"
    h0::Float64 = 5960
    "radius of conical mountain"
    a::Float64 = 20.0
    "center of mountain long coord, shifted by 180 compared to the paper, 
    because our λ ∈ [-180, 180] (in the paper it was 270, with λ ∈ [0, 360])"
    λc::Float64 = 90.0
    "latitude coordinate for center of mountain"
    ϕc::Float64 = 30.0
    "mountain peak height"
    h_s0::Float64 = 2e3
end

function surface_topography(space::Spaces.SpectralElementSpace2D, test::MountainTest)
    FT = Spaces.undertype(space)
    a = FT(test.a)
    λc = FT(test.λc)
    ϕc = FT(test.ϕc)
    h_s0 = FT(test.h_s0)
    map(Fields.coordinate_field(space)) do coord
        ϕ = coord.lat
        λ = coord.long
        r = sqrt(min(a^2, (λ - λc)^2 + (ϕ - ϕc)^2)) # positive branch
        h_s0 * (1 - r / a)
    end
end
function initial_condition(
    space,
    test::MountainTest,
)
    # steady-state and mountain test cases share the same initial condition
    initial_condition(space, SteadyStateTest(test.params, test.u0, test.h0))
end


"""
    RossbyHaurwitzTest <: AbstractSphereTestCase

A fourth one, called "rossby_haurwitz", reproduces Test Case 6 in the same
reference paper. It represents the solution of the nonlinear barotropic
vorticity equation on the sphere

https://doi.org/10.1016/S0021-9991(05)80016-6

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct RossbyHaurwitzTest <: AbstractSphereTestCase
    "Physical parameters"
    params::SphericalParameters = SphericalParameters()
    "velocity amplitude parameter"
    a::Float64 = 4.0
    "peak of analytic height field"
    h0::Float64 = 8.0e3
    "vorticity amplitude parameter (1/sec)"
    ω::Float64 = 7.848e-6
    "vorticity amplitude parameter (1/sec)"
    K::Float64 = 7.848e-6
end
RossbyHaurwitzTest(α::FT) where {FT} =
    RossbyHaurwitzTest{FT, PhysicalParameters{FT}}(; α = α)
    
function initial_condition(space, test::RossbyHaurwitzTest)
    FT = Spaces.undertype(space)
    a = FT(test.a)
    h0 = FT(test.h0)
    ω = FT(test.ω)
    K = FT(test.K)
    R = FT(test.params.R)
    Ω = FT(test.params.Ω)
    g = FT(test.params.g)

    Y = map(Fields.local_geometry_field(space)) do local_geometry
        coord = local_geometry.coordinates
        ϕ = coord.lat
        λ = coord.long

        A =
            ω / 2 * (2 * Ω + ω) * cosd(ϕ)^2 +
            1 / 4 *
            K^2 *
            cosd(ϕ)^(2 * a) *
            ((a + 1) * cosd(ϕ)^2 + (2 * a^2 - a - 2) - 2 * a^2 * cosd(ϕ)^-2)
        B =
            2 * (Ω + ω) * K / (a + 1) / (a + 2) *
            cosd(ϕ)^a *
            ((a^2 + 2 * a + 2) - (a + 1)^2 * cosd(ϕ)^2)
        C = 1 / 4 * K^2 * cosd(ϕ)^(2 * a) * ((a + 1) * cosd(ϕ)^2 - (a + 2))

        h =
            h0 +
            (R^2 * A + R^2 * B * cosd(a * λ) + R^2 * C * cosd(2 * a * λ)) / g

        uλ =
            R * ω * cosd(ϕ) +
            R * K * cosd(ϕ)^(a - 1) * (a * sind(ϕ)^2 - cosd(ϕ)^2) * cosd(a * λ)
        uϕ = -R * K * a * cosd(ϕ)^(a - 1) * sind(ϕ) * sind(a * λ)


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
    BarotropicInstabilityTest{FT, P} <: AbstractSphereTestCase

A fifth one, called "barotropic_instability", reproduces the test case in
Galewsky et al, "An initial-value problem for testing numerical models of
the global shallow-water equations", 2004 (also in Sec. 7.6 of Ullirch et al,
"High-order ﬁnite-volume methods for the shallow-water equations on
the sphere", 2010). This test case consists of a zonal jet with compact
support at a latitude of 45°. A small height disturbance is then added,
which causes the jet to become unstable and collapse into a highly vortical
structure.

https://doi.org/10.3402/tellusa.v56i5.14436
https://doi.org/10.1016/j.jcp.2010.04.044

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct BarotropicInstabilityTest{FT, P} <: AbstractSphereTestCase
    "Physical parameters"
    params::SphericalParameters = SphericalParameters()
    "maximum zonal velocity"
    u_max::Float64 = 80.0
    "mountain shape parameters"
    αₚ::Float64 = 19.09859
    "mountain shape parameters"
    βₚ::Float64 = 3.81971
    "peak of balanced height field from Tempest 
    https://github.com/paullric/tempestmodel/blob/master/test/shallowwater_sphere/BarotropicInstabilityTest.cpp#L86"
    h0::Float64 = 10158.18617
    "local perturbation peak height"
    h_hat::Float64 = 120.0
    "southern jet boundary"
    ϕ₀::Float64 = 25.71428
    "northern jet boundary"
    ϕ₁::Float64 = 64.28571
    "height perturbation peak location"
    ϕ₂::Float64 = 45.0
    "zonal velocity decay parameter"
    eₙ::Float64 = exp(-4.0 / (deg2rad(ϕ₁) - deg2rad(ϕ₀))^2)
end

function initial_condition(space, test::BarotropicInstabilityTest)
    FT = Spaces.undertype(space)
    u_max = FT(test.u_max)
    αₚ = FT(test.αₚ)
    βₚ = FT(test.βₚ)
    h0 = FT(test.h0)
    h_hat = FT(test.h_hat)
    ϕ₀ = FT(test.ϕ₀)
    ϕ₁ = FT(test.ϕ₁)
    ϕ₂ = FT(test.ϕ₂)
    eₙ = FT(test.eₙ)
    R = FT(test.params.R)
    Ω = FT(test.params.Ω)
    g = FT(test.params.g)
    α = FT(test.params.α)

    Y = map(Fields.local_geometry_field(space)) do local_geometry
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
        h =
            h0 - (R / g) * (pi / 180.0) * QuadGK.quadgk(h_int, -90.0, ϕprime)[1]

        if λ > 0.0
            λ -= 360.0
        end
        if λ < -360.0 || λ > 0.0
            @info "Invalid longitude value"
        end

        # Add height perturbation
        h += h_hat * cosd(ϕ) * exp(-(λ^2 / αₚ^2) - ((ϕ₂ - ϕ)^2 / βₚ^2))

        uλ = uλprime(ϕprime)
        uϕ = uϕprime

        u = Geometry.transform(
            Geometry.Covariant12Axis(),
            Geometry.UVVector(uλ, uϕ),
            local_geometry,
        )

        return (h = h, u = u)
    end
    return Y
end