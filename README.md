# ClimaShallowWater.jl

[![Build status](https://badge.buildkite.com/4c3b1db756d1987ab7920a17579c37142623d356b7b77fe670.svg?branch=main)](https://buildkite.com/clima/climashallowwater-ci)

This is an implementation of the [shallow water equations](https://en.wikipedia.org/wiki/Shallow_water_equations) using [ClimaCore.jl](https://github.com/CliMA/ClimaCore.jl).

## Installation

You will need an installation of Julia (https://julialang.org/downloads/), and then run
```sh
julia --project -e 'using Pkg; Pkg.instantiate()'
```
to download and install any necessary dependencies.

## Running

The main driver file is `shallowwater`. Run
```sh
./shallowwater --help
```
to see the arguments.