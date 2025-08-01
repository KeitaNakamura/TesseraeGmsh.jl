using TesseraeGmsh
using Tesserae
using Test

@testset "readmsh" begin
    meshes = TesseraeGmsh.readmsh("sphere.msh")
    # body
    mesh_body = meshes["body"]
    particles = generate_particles(@NamedTuple{x::Vec{3,Float64}, detJdΩ::Float64}, mesh_body)
    weights = generate_interpolation_weights(mesh_body, size(particles))
    feupdate!(weights, mesh_body; volume = particles.detJdΩ)
    @test sum(particles.detJdΩ) ≈ 4π/3 rtol=0.01
    # surface
    mesh_surface = meshes["surface"]
    particles = generate_particles(@NamedTuple{x::Vec{3,Float64}, detJdΩ::Float64}, mesh_surface)
    weights = generate_interpolation_weights(mesh_surface, size(particles))
    feupdate!(weights, mesh_surface; area = particles.detJdΩ)
    @test sum(particles.detJdΩ) ≈ 4π rtol=0.001
end
