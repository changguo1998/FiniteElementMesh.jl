
export triangle_norm_out, tetra_norm_out

function _triangle_norm_out(p1::Vector, p2::Vector, pref::Vector)
    r12 = p2 - p1
    r31 = p1 - pref
    n = normalize([r12[2], -r12[1]])
    if dot(r31, n) < 0.0
        n .*= -1.0
    end
    return n
end

"""
```julia
triangle_norm_out(mesh::Mesh2D [, aux::Mesh2Daux]) -> Mesh2Ddata(:normvec, :edge)
```
"""
function triangle_norm_out(mesh::Mesh2D, aux::Union{Nothing,Mesh2Daux}=nothing)
    if isnothing(aux)
        aux = auxillary_relation(mesh)
    end
    buffer = zeros(2, 3, nface(mesh))
    for f = axes(mesh.face2edge, 2)
        facevert_list = aux.face2vert[:, f]
        for el = axes(mesh.face2edge, 1)
            edgevert_list = mesh.edge2vert[:, mesh.face2edge[el, f]]
            p1 = mesh.vertex[:, edgevert_list[1]]
            p2 = mesh.vertex[:, edgevert_list[2]]
            p3 = mesh.vertex[:, first(setdiff(facevert_list, edgevert_list))]
            buffer[:, el, f] = _triangle_norm_out(p1, p2, p3)
        end
    end
    return Mesh2Ddata((:normvec, :edge, :face), buffer)
end

function _tetra_norm_out(p1::Vector, p2::Vector, p3::Vector, pref::Vector)
    r12 = p2 - p1
    r13 = p3 - p1
    r41 = p1 - pref
    n = normalize(cross(r12, r13))
    if dot(n, r41) < 0.0
        n .*= -1.0
    end
    return n
end

"""
```julia
tetra_norm_out(mesh::Mesh3D [, aux::Mesh3Daux]) -> Mesh3Ddata(:normvec, :face, :cell)
```
"""
function tetra_norm_out(mesh::Mesh3D, aux::Union{Nothing,Mesh3Daux}=nothing)
    if isnothing(aux)
        aux = auxillary_relation(mesh)
    end

    buffer = zeros(3, 4, ncell(mesh))
    for c = axes(mesh.cell2face, 2)
        cellvert_list = aux.cell2vert[:, c]
        for fl = axes(mesh.cell2face, 1)
            facevert_list = aux.face2vert[:, mesh.cell2face[fl, c]]
            p1 = mesh.vertex[:, facevert_list[1]]
            p2 = mesh.vertex[:, facevert_list[2]]
            p3 = mesh.vertex[:, facevert_list[3]]
            p4 = mesh.vertex[:, first(setdiff(cellvert_list, facevert_list))]
            buffer[:, fl, c] = _tetra_norm_out(p1, p2, p3, p4)
        end
    end
    return Mesh3Ddata((:normv, :face, :cell), buffer)
end
