# TesseraeGmsh.jl

## Installation

```julia
pkg> add https://github.com/KeitaNakamura/TesseraeGmsh.jl.git
```

## Usage

```julia
julia> using TesseraeGmsh

julia> meshes = TesseraeGmsh.readmsh("test/sphere.msh")
Info    : Reading 'test/sphere.msh'...
Info    : 7 entities
Info    : 10888 nodes
Info    : 64098 elements
Info    : Done reading 'test/sphere.msh'
Dict{String, Tesserae.UnstructuredMesh} with 2 entries:
  "body"    => Vec{3, Float64}[[6.12323e-17, -1.49976e-32, 1.0], …
  "surface" => Vec{3, Float64}[[6.12323e-17, -1.49976e-32, 1.0], …
```
