
function _read_arr(io::IO, T::Type)
    nd = read(io, Int)
    sz = zeros(Int, nd)
    read!(io, sz)
    arr = zeros(T, sz)
    read!(io, arr)
    return arr
end

function _write_arr(io::IO, arr::Array)
    if eltype(arr) <: Int
        write(io, ndims(arr))
        write(io, Int.(size(arr)))
        write(io, Int.(arr))
    elseif eltype(arr) <: Float64
        write(io, ndims(arr))
        write(io, Int.(size(arr)))
        write(io, Float64.(arr))
    end
end

export FEMesh, FEMeshData, show, dimension, nvertex, nedge, nface, ncell

abstract type FEMesh <: Any end
abstract type FEMeshData <: Any end

#
# Mesh1D
#

export Mesh1D, Mesh1Daux, Mesh1Ddata

struct Mesh1D <: FEMesh
    vertex::Vector{Float64}
    edge2vert::Matrix{Int}
end

dimension(mesh::Mesh1D) = 1
nvertex(mesh::Mesh1D) = length(mesh.vertex)
nedge(mesh::Mesh1D) = size(mesh.edge2vert, 2)

function show(io::IO, ::MIME"text/plain", m::Mesh1D)
    println(io, "vertex: ", size(m.vertex))
    print(io, "edge2vert: ", size(m.edge2vert))
end

function write(io::IO, mesh::Mesh1D)
    _write_arr(io, mesh.vertex)
    _write_arr(io, mesh.edge2vert)
end

function read(io::IO, ::Type{Mesh1D})
    v = _read_arr(io, Float64)
    d2v = _read_arr(io, Int)
    return Mesh1D(vec(v), d2v)
end

struct Mesh1Daux
    # type B: low dim -> high dim
    vert2edge::Matrix{Int}

    # type C: low dim -> high dim id
    vert2edge_l::Matrix{Int}
end

function show(io::IO, ::MIME"text/plain", m::Mesh1Daux)
    print(io, "vert2edge: ", size(m.vert2edge))
end

struct Mesh1Ddata <: FEMeshData
    dimName::NTuple
    data::Array
end

function show(io::IO, ::MIME"text/plain", m::Mesh1Ddata)
    s = size(m.data)
    for i = eachindex(s)
        println(io, m.dimName[i], " ", s[i])
    end
end

#
# Mesh2D
#

export Mesh2D, Mesh2Daux, Mesh2Ddata

struct Mesh2D <: FEMesh
    vertex::Matrix{Float64}
    edge2vert::Matrix{Int}
    face2edge::Matrix{Int}
end

dimension(mesh::Mesh2D) = 2
nvertex(mesh::Mesh2D) = size(mesh.vertex, 2)
nedge(mesh::Mesh2D) = size(mesh.edge2vert, 2)
nface(mesh::Mesh2D) = size(mesh.face2edge, 2)

function show(io::IO, ::MIME"text/plain", m::Mesh2D)
    println(io, "vertex: ", size(m.vertex))
    println(io, "edge2vert: ", size(m.edge2vert))
    print(io, "face2edge: ", size(m.face2edge))
end

function write(io::IO, mesh::Mesh2D)
    _write_arr(io, mesh.vertex)
    _write_arr(io, mesh.edge2vert)
    _write_arr(io, mesh.face2edge)
end

function read(io::IO, ::Type{Mesh2D})
    v = _read_arr(io, Float64)
    d2v = _read_arr(io, Int)
    f2d = _read_arr(io, Int)
    return Mesh2D(v, d2v, f2d)
end

struct Mesh2Daux
    # type A: high dim -> low dim
    face2vert::Matrix{Int}

    # type B: low dim -> high dim
    vert2edge::Matrix{Int}
    vert2face::Matrix{Int}
    edge2face::Matrix{Int}

    # type C: low dim -> high dim local id
    vert2edge_l::Matrix{Int}
    vert2face_l::Matrix{Int}
    edge2face_l::Matrix{Int}
end

function show(io::IO, ::MIME"text/plain", m::Mesh2Daux)
    println(io, "face2vert: ", size(m.face2vert))
    println(io, "edge2face: ", size(m.edge2face))
    println(io, "vert2edge: ", size(m.vert2edge))
    print(io, "vert2face: ", size(m.vert2face))
end

struct Mesh2Ddata <: FEMeshData
    dimName::NTuple
    data::Array
end

function show(io::IO, ::MIME"text/plain", m::Mesh2Ddata)
    s = size(m.data)
    for i = eachindex(s)
        println(io, m.dimName[i], " ", s[i])
    end
end

#
# Mesh3D
#

export Mesh3D, Mesh3Daux, Mesh3Ddata

struct Mesh3D <: FEMesh
    vertex::Matrix{Float64}
    edge2vert::Matrix{Int}
    face2edge::Matrix{Int}
    cell2face::Matrix{Int}
end

dimension(mesh::Mesh3D) = 3
nvertex(mesh::Mesh3D) = size(mesh.vertex, 2)
nedge(mesh::Mesh3D) = size(mesh.edge2vert, 2)
nface(mesh::Mesh3D) = size(mesh.face2edge, 2)
ncell(mesh::Mesh3D) = size(mesh.cell2face, 2)

function show(io::IO, ::MIME"text/plain", m::Mesh3D)
    println(io, "vertex: ", size(m.vertex))
    println(io, "edge2vert: ", size(m.edge2vert))
    println(io, "face2edge: ", size(m.face2edge))
    print(io, "cell2face: ", size(m.cell2face))
end

function write(io::IO, mesh::Mesh3D)
    _write_arr(io, mesh.vertex)
    _write_arr(io, mesh.edge2vert)
    _write_arr(io, mesh.face2edge)
    _write_arr(io, mesh.cell2face)
end

function read(io::IO, ::Type{Mesh3D})
    v = _read_arr(io, Float64)
    d2v = _read_arr(io, Int)
    f2d = _read_arr(io, Int)
    e2f = _read_arr(io, Int)
    return Mesh3D(v, d2v, f2d, e2f)
end

struct Mesh3Daux
    # type A: high dim -> low dim
    face2vert::Matrix{Int}
    cell2edge::Matrix{Int}
    cell2vert::Matrix{Int}

    # type B: low dim -> high dim
    vert2edge::Matrix{Int}
    vert2face::Matrix{Int}
    vert2cell::Matrix{Int}
    edge2face::Matrix{Int}
    edge2cell::Matrix{Int}
    face2cell::Matrix{Int}

    # type C: low dim -> high dim local id
    vert2edge_l::Matrix{Int}
    vert2face_l::Matrix{Int}
    vert2cell_l::Matrix{Int}
    edge2face_l::Matrix{Int}
    edge2cell_l::Matrix{Int}
    face2cell_l::Matrix{Int}
end

function show(io::IO, ::MIME"text/plain", m::Mesh3Daux)
    println(io, "face2vert: ", size(m.face2vert))
    println(io, "cell2edge: ", size(m.cell2edge))
    println(io, "cell2vert: ", size(m.cell2vert))
    println(io, "vert2edge: ", size(m.vert2edge))
    println(io, "vert2face: ", size(m.vert2face))
    println(io, "vert2cell: ", size(m.vert2cell))
    println(io, "edge2face: ", size(m.edge2face))
    println(io, "edge2cell: ", size(m.edge2cell))
    print(io, "face2cell: ", size(m.face2cell))
end

struct Mesh3Ddata <: FEMeshData
    dimName::NTuple
    data::Array
end

function show(io::IO, ::MIME"text/plain", m::Mesh3Ddata)
    s = size(m.data)
    for i = eachindex(s)
        println(io, m.dimName[i], " ", s[i])
    end
end
