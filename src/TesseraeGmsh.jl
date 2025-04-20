module TesseraeGmsh

using Tesserae
using StaticArrays

include("gmsh.jl")
import .GmshReader

function from_gmsh(nodeset::GmshReader.NodeSet)
    dim = nodeset.dim
    map(v -> Vec{dim}(i -> v[i]), nodeset.coord)
end

function from_gmsh_shape(shape::String)
    shape == "Line 2" && return Tesserae.Line2()
    shape == "Line 3" && return Tesserae.Line3()
    shape == "Triangle 3" && return Tesserae.Tri3()
    shape == "Triangle 6" && return Tesserae.Tri6()
    shape == "Tetrahedron 4" && return Tesserae.Tet4()
    shape == "Tetrahedron 10" && return Tesserae.Tet10()
    shape == "Quadrilateral 4" && return Tesserae.Quad4()
    shape == "Quadrilateral 9" && return Tesserae.Quad9()
    shape == "Hexahedron 8" && return Tesserae.Hex8()
    shape == "Hexahedron 27" && return Tesserae.Hex27()
    error("\"$shape\" is not supported yet")
end

from_gmsh_connectivity(::Tesserae.Line2) = SVector(1,2)
from_gmsh_connectivity(::Tesserae.Line3) = SVector(1,2,3)
from_gmsh_connectivity(::Tesserae.Tri3)  = SVector(1,2,3)
from_gmsh_connectivity(::Tesserae.Tri6)  = SVector(1,2,3,4,6,5)
from_gmsh_connectivity(::Tesserae.Tet4)  = SVector(1,2,3,4)
from_gmsh_connectivity(::Tesserae.Tet10) = SVector(1,2,3,4,5,7,8,6,9,10)
from_gmsh_connectivity(::Tesserae.Quad4) = SVector(1,2,3,4)
from_gmsh_connectivity(::Tesserae.Quad9) = SVector(1,2,3,4,5,6,7,8,9)
from_gmsh_connectivity(::Tesserae.Hex8)  = SVector(1,2,3,4,5,6,7,8)
from_gmsh_connectivity(::Tesserae.Hex27) = SVector(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27)

function from_gmsh(elementset::GmshReader.ElementSet)
    shape = from_gmsh_shape(elementset.elementname)
    connectivities = map(elementset.connectivities) do conn
        @assert Tesserae.nlocalnodes(shape) == length(conn)
        conn[from_gmsh_connectivity(shape)]
    end
    shape, connectivities
end

function from_gmsh(phygroup::GmshReader.PhysicalGroup)
    # elements
    dict = Dict{Tesserae.Shape, Vector{<: SVector}}()
    for entitiy in phygroup.entities
        for elementset in entitiy
            shape, connectivities = from_gmsh(elementset)
            append!(get!(dict, shape, SVector{Tesserae.nlocalnodes(shape), Int}[]), connectivities)
        end
    end
    only(dict) # support only unique shape
end

function from_gmsh(gmshfile::GmshReader.GmshFile)
    dim = gmshfile.nodeset.dim
    # nodes
    nodes = from_gmsh(gmshfile.nodeset)
    # physicalgroups
    Dict{String, UnstructuredMesh}(Iterators.map(gmshfile.physicalgroups) do (name, phygroup)
        shape, connectivities = from_gmsh(phygroup)
        name => UnstructuredMesh(shape, nodes, connectivities)
    end)
end

function readmsh(filename::String)
    file = GmshReader.readmsh(filename; outward_surface_normals = true)
    from_gmsh(file)
end

end # module TesseraeGmsh
