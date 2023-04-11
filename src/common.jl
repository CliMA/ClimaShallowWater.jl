abstract type AbstractSphereTestCase end

function create_space(
    context,
    test::AbstractSphereTestCase;
    float_type = Float64,
    panel_size = 9,
    poly_nodes = 4,
)
    domain = Domains.SphereDomain(float_type(test.params.R))
    mesh = Meshes.EquiangularCubedSphere(domain, panel_size)
    quad = Spaces.Quadratures.GLL{poly_nodes}()
    topology = Topologies.Topology2D(context, mesh)
    space = Spaces.SpectralElementSpace2D(topology, quad)
    return space
end

function approx_resolution(space::Spaces.SpectralElementSpace2D)
    FT = Spaces.undertype(space)
    mesh = space.topology.mesh
    mesh_scale = FT(sqrt(4 * pi / 6)) * mesh.domain.radius / mesh.ne

    quad = Spaces.quadrature_style(space)
    Nq = Spaces.Quadratures.degrees_of_freedom(quad)
    return mesh_scale / (Nq - 1)
end


"""
    surface_height_field(space, test::AbstractSphereTestCase)

Construct a field on `space` of the surface topography height `h_s` for the
test case `test`.

Set to be a field of zeros if not otherwise defined.
"""
function surface_height_field(
    space::Spaces.SpectralElementSpace2D,
    test::AbstractSphereTestCase,
)
    return zeros(space)
end

"""
    coriolis_field(space, test::AbstractSphereTestCase)

Construct a field on `space` of the angular coriolis velocity for the test case `test`.

By default it is
```math
f = 2Ω (\\sin(\\phi) \\cos(\\alpha) - \\cos(\\lambda) \\cos(\\phi) \\sin(\\alpha)
```
where ``\\alpha`` is the angle between the north pole and the center of the top, and
``\\lambda`` and ``\\phi`` are the longitude and latitude of the grid point.
"""
function coriolis_field(
    space::Spaces.SpectralElementSpace2D,
    test::AbstractSphereTestCase,
)
    coordinates = Fields.coordinate_field(space)
    FT = Spaces.undertype(space)
    ϕ = coordinates.lat
    λ = coordinates.long

    Ω = FT(test.params.Ω)
    α = FT(test.params.α)

    f = @. Geometry.Contravariant3Vector(
        Geometry.WVector(
            2 * Ω * (sind(ϕ) * cosd(α) - cosd(λ) * cosd(ϕ) * sind(α)),
        ),
    )
    return f
end

function hyperdiffusion_coefficient(
    space::Spaces.SpectralElementSpace2D,
    test::AbstractSphereTestCase,
)
    FT = Spaces.undertype(space)
    ν = test.params.ν
    c = sqrt(test.params.g * test.h0)
    D₄ = ν * c * approx_resolution(space)^3
    return FT(D₄)
end

function auxiliary_state(Y, test::AbstractSphereTestCase)
    h = Y.h
    u = Y.u
    space = axes(h)

    FT = Spaces.undertype(space)

    f = coriolis_field(space, test)
    h_s = surface_height_field(space, test)

    h_buffer = Spaces.create_dss_buffer(h)
    u_buffer = Spaces.create_dss_buffer(u)

    D₄ = hyperdiffusion_coefficient(space, test)
    (; g = FT(test.params.g), D₄, f, h_s, h_buffer, u_buffer)
end

function dss!(Y, p)
    (; h, u) = Y
    (; h_buffer, u_buffer) = p
    NVTX.@range "dss" begin
        Spaces.weighted_dss_start!(h, h_buffer)
        Spaces.weighted_dss_start!(u, u_buffer)

        Spaces.weighted_dss_internal!(h, h_buffer)
        Spaces.weighted_dss_internal!(u, u_buffer)

        Spaces.weighted_dss_ghost!(h, h_buffer)
        Spaces.weighted_dss_ghost!(u, u_buffer)
    end
end

"""
    initial_condition(space, test::AbstractSphereTestCase)

Construct the initial state for the test case `test` on the `space` provided.
"""
function initial_condition(
    space::Spaces.SpectralElementSpace2D,
    test::AbstractSphereTestCase,
)
    h = initial_height(space, test)
    u = initial_velocity(space, test)
    return Fields.FieldVector(h = h, u = u)
end



function tendency!(dYdt, y, p, t)
    NVTX.@range "rhs!" begin
        (; f, h_s, D₄, g) = p

        div = Operators.Divergence()
        wdiv = Operators.WeakDivergence()
        grad = Operators.Gradient()
        wgrad = Operators.WeakGradient()
        curl = Operators.Curl()
        wcurl = Operators.WeakCurl()

        # Compute Laplacians of h and u
        NVTX.@range "hyperdiffusion" begin
            @. dYdt.h = wdiv(grad(y.h))
            @. dYdt.u =
                wgrad(div(y.u)) - Geometry.Covariant12Vector(
                    wcurl(Geometry.Covariant3Vector(curl(y.u))),
                )

            dss!(dYdt, p)
        end

        NVTX.@range "tendency" begin
            @. dYdt.h = -D₄ * wdiv(grad(dYdt.h))
            # split to avoid "device kernel image is invalid (code 200, ERROR_INVALID_IMAGE)"
            @. dYdt.u =
                wgrad(div(dYdt.u)) - Geometry.Covariant12Vector(
                    wcurl(Geometry.Covariant3Vector(curl(dYdt.u))),
                )
            @. dYdt.u = -D₄ * dYdt.u
            @. begin
                dYdt.h += -wdiv(y.h * y.u)
                dYdt.u += -grad(g * (y.h + h_s) + norm(y.u)^2 / 2) #+
                dYdt.u += y.u × (f + curl(y.u))
            end
            dss!(dYdt, p)
        end
    end
    return dYdt
end
