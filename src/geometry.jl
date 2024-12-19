
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

function triangle_norm_out(mesh::Mesh2D, mesha::Union{Nothing,Mesh2Daux}=nothing)
    if isnothing(mesha)
        mesha = auxillary_relation(mesh)
    end
    Nface = nface(mesh)
    buffer = zeros(2, 3, Nface)
    for iface = axes(mesh.face2edge, 2)
        verts = mesha.face2vert[:, iface]
        for idl = axes(mesh.face2edge, 1)
            edgeverts = mesh.edge2vert[:, mesh.face2edge[idl, iface]]
            p1 = mesh.vertex[:, edgeverts[1]]
            p2 = mesh.vertex[:, edgeverts[2]]
            ip3 = setdiff(verts, edgeverts)[1]
            p3 = mesh.vertex[:, ip3]
            buffer[:, idl, iface] = _triangle_norm_out(p1, p2, p3)
        end
    end
    return Mesh2Ddata((:normv, :edge, :face), buffer)
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

function tetra_norm_out(mesh::Mesh3D, mesha::Union{Nothing,Mesh3Daux}=nothing)
    if isnothing(mesha)
        mesha = auxillary_relation(mesh)
    end

    Ncell = ncell(mesh)
    buffer = zeros(3, 4, Ncell)
    for icell = axes(mesh.cell2face, 2)
        vertexs = mesha.cell2vert[:, icell]
        for ifl = axes(mesh.cell2face, 1)
            vertexf = mesha.face2vert[:, mesh.cell2face[ifl, icell]]
            p1 = mesh.vertex[:, vertexf[1]]
            p2 = mesh.vertex[:, vertexf[2]]
            p3 = mesh.vertex[:, vertexf[3]]
            ip4 = first(setdiff(vertexs, vertexf))
            p4 = mesh.vertex[:, ip4]
            buffer[:, ifl, icell] = _tetra_norm_out(p1, p2, p3, p4)
        end
    end
    return Mesh3Ddata((:normv, :face, :cell), buffer)
end