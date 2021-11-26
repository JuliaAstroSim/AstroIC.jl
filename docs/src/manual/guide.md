# Package Guide

## Installation

From the Julia REPL, type `]` to enter the Pkg REPL mode and run
```julia
pkg> add AstroIC
```
or add from git repository
```julia
pkg> add https://github.com/JuliaAstroSim/AstroIC.jl
```

Test the package by
```julia
pkg> test AstroIC
```

## Basic Usage

```@repl guide
using AstroIC

config = PlummerStarCluster()
data = generate(config)
```