# ClimaShallowWater.jl

ClimaShallowWater.jl is an implementation of the [shallow water equations](https://en.wikipedia.org/wiki/Shallow_water_equations) using [ClimaCore.jl](https://github.com/CliMA/ClimaCore.jl).

```@meta
CurrentModule = ClimaShallowWater
```

## Equations

This solves the shallow-water equations in vector-invariant form: there are two prognostic variables _height_ from the surface $h$ and _velocity_ $\boldsymbol{u}$, with tendency terms

```math
\frac{\partial h}{\partial t} = - \nabla \cdot (h \boldsymbol{u})
```

```math
\frac{\partial \boldsymbol{u}}{\partial t} = - \nabla \left[ g (h + h_s) + \tfrac{1}{2} \|\boldsymbol{u}\|^2 \right] + \boldsymbol{u} \times [\boldsymbol{f} + \nabla \times (\boldsymbol{u})]
```
where $h_s$ is the surface height, and $\boldsymbol{f}$ is the Coriolis parameter.
