# To make differences of same data between local(new) mesh and global(old) mesh,
# the variable is appended with suffix `_l` and `_g`
# suffix `_l2g` and `_g2l` is also added to show the usage for index converting

_partition_match_vec(vec1::Vector{Int}, vec2::Vector{Int}) = vec1 == vec2
_partition_match_vec(vec1::Vector{Float64}, vec2::Vector{Float64}) = vec1 == vec2

function _partition_match_vec_ignore_range(vec1::Vector{Int}, vec2::Vector{Int})
    v1 = sort(vec1)
    v2 = sort(vec2)
    return v1 == v2
end

function _partition_match_col(mf::Function, relation::Matrix{Int}, reference::Vector{Int})
    for ic in axes(relation, 2)
        if mf(relation[:, ic], reference)
            return ic
        end
    end
    return 0
end

function _partition_index_array(arr::Array{<:Any}, idx::Tuple)
    idxvec = map(idx) do id
        if typeof(id) <: Integer
            Int.(collect(1:id))
        elseif typeof(id) <: AbstractVector{<:Integer}
            Int.(collect(id))
        else
            return Int[]
        end
    end
    newdim = length(idxvec)
    newsize = Tuple(length.(idxvec))
    newdata = similar(arr, newsize)
    for lidx = CartesianIndices(newsize)
        gidx = zeros(Int, newdim)
        for idim = eachindex(newsize)
            gidx[idim] = idxvec[idim][lidx[idim]]
        end
        newdata[lidx] = arr[CartesianIndex(Tuple(gidx))]
    end
    return newdata
end

function _partition_index_l2g_to_g2l(index_l2g::Vector{Int}, n_global::Int)
    index_g2l = zeros(Int, n_global)
    foreach(eachindex(index_l2g)) do il
        ig = index_l2g[il]
        index_g2l[ig] = il
    end
    return index_g2l
end

function _partition_get_Xl2g_from_Yl2g(Y2X::Matrix{Int}, Ysub_l2g::Vector{Int}, nx_g::Int)
    flagX_g = falses(nx_g)
    for iy_g = Ysub_l2g
        for ix_g = Y2X[:, iy_g]
            flagX_g[ix_g] = true
        end
    end
    return findall(flagX_g)
end

function _partition_get_Yl2g_from_Xl2g(Y2X::Matrix{Int}, Xsub_l2g::Vector{Int})
    flagmatch = map(ix -> ix in Xsub_l2g, Y2X)
    flagY = all(flagmatch, dims=1) |> vec
    return findall(flagY)
end

function _partition_get_local(Y2X_g::Matrix{Int}, X_g2l::Vector{Int}, Y_l2g::Vector{Int})
    return X_g2l[Y2X_g[:, Y_l2g]]
end

function _partition_get_global(Y2X_l::Matrix{Int}, X_l2g::Vector{Int})
    return X_l2g[Y2X_l]
end

#
#
#

export submesh

function submesh_v(mesh::Mesh1D, vertex_l2g::AbstractVector{<:Integer})
    s_vertices = mesh.vertex[vertex_l2g]

    vertex_g2l = _partition_index_l2g_to_g2l(vertex_l2g, nvertex(mesh))
    edge_l2g = _partition_get_Yl2g_from_Xl2g(mesh.edge2vert, vertex_l2g)
    s_edge2vert = _partition_get_local(mesh.edge2vert, vertex_g2l, edge_l2g)
    return Mesh1D(s_vertices, s_edge2vert)
end

function submesh_e(mesh::Mesh1D, edge_l2g::AbstractVector{<:Integer})
    vertex_l2g = _partition_get_Xl2g_from_Yl2g(mesh.edge2vert, edge_l2g, nvertex(mesh))
    vertex_g2l = _partition_index_l2g_to_g2l(vertex_l2g, nvertex(mesh))

    s_vertices = mesh.vertex[vertex_l2g]
    s_edge2vert = _partition_get_local(mesh.edge2vert, vertex_g2l, edge_l2g)
    return Mesh1D(s_vertices, s_edge2vert)
end

function submesh(mesh::Mesh1D, var::Symbol, l2g::AbstractVector{<:Integer})
    if var == :vertex
        return submesh_v(mesh, l2g)
    elseif var == :edge
        return submesh_e(mesh, l2g)
    else
        error("var: $(var) not supported, availables are: :vertex, :edge")
    end
end

function submesh_v(mesh::Mesh2D, vertex_l2g::AbstractVector{<:Integer})
    s_vertices = mesh.vertex[:, vertex_l2g]

    vertex_g2l = _partition_index_l2g_to_g2l(vertex_l2g, nvertex(mesh))
    edge_l2g = _partition_get_Yl2g_from_Xl2g(mesh.edge2vert, vertex_l2g)
    s_edge2vert = _partition_get_local(mesh.edge2vert, vertex_g2l, edge_l2g)

    edge_g2l = _partition_index_l2g_to_g2l(edge_l2g, nedge(mesh))
    face_l2g = _partition_get_Yl2g_from_Xl2g(mesh.face2edge, edge_l2g)
    s_face2edge = _partition_get_local(mesh.face2edge, edge_g2l, face_l2g)
    return Mesh2D(s_vertices, s_edge2vert, s_face2edge)
end

function submesh_e(mesh::Mesh2D, edge_l2g::AbstractVector{<:Integer})
    face_l2g = _partition_get_Yl2g_from_Xl2g(mesh.face2edge, edge_l2g)
    edge_g2l = _partition_index_l2g_to_g2l(edge_l2g, nedge(mesh))
    s_face2edge = _partition_get_local(mesh.face2edge, edge_g2l, face_l2g)

    vertex_l2g = _partition_get_Xl2g_from_Yl2g(mesh.edge2vert, edge_l2g, nvertex(mesh))
    vertex_g2l = _partition_index_l2g_to_g2l(vertex_l2g, nvertex(mesh))
    s_edge2vert = _partition_get_local(mesh.edge2vert, vertex_g2l, edge_l2g)

    s_vertices = mesh.vertex[:, vertex_l2g]
    return Mesh2D(s_vertices, s_edge2vert, s_face2edge)
end

function submesh_f(mesh::Mesh2D, face_l2g::AbstractVector{<:Integer})
    edge_l2g = _partition_get_Xl2g_from_Yl2g(mesh.face2edge, face_l2g, nedge(mesh))
    edge_g2l = _partition_index_l2g_to_g2l(edge_l2g, nedge(mesh))
    s_face2edge = _partition_get_local(mesh.face2edge, edge_g2l, face_l2g)

    vertex_l2g = _partition_get_Xl2g_from_Yl2g(mesh.edge2vert, edge_l2g, nvertex(mesh))
    vertex_g2l = _partition_index_l2g_to_g2l(vertex_l2g, nvertex(mesh))
    s_edge2vert = _partition_get_local(mesh.edge2vert, vertex_g2l, edge_l2g)

    s_vertices = mesh.vertex[:, vertex_l2g]
    return Mesh2D(s_vertices, s_edge2vert, s_face2edge)
end

function submesh(mesh::Mesh2D, var::Symbol, l2g::AbstractVector{<:Integer})
    if var == :vertex
        return submesh_v(mesh, l2g)
    elseif var == :edge
        return submesh_e(mesh, l2g)
    elseif var == :face
        return submesh_f(mesh, l2g)
    else
        error("var: $(var) not supported, availables are: :vertex, :edge, :face")
    end
end

function submesh_v(mesh::Mesh3D, vertex_l2g::AbstractVector{<:Integer})
    s_vertices = mesh.vertex[:, vertex_l2g]

    vertex_g2l = _partition_index_l2g_to_g2l(vertex_l2g, nvertex(mesh))
    edge_l2g = _partition_get_Yl2g_from_Xl2g(mesh.edge2vert, vertex_l2g)
    s_edge2vert = _partition_get_local(mesh.edge2vert, vertex_g2l, edge_l2g)

    edge_g2l = _partition_index_l2g_to_g2l(edge_l2g, nedge(mesh))
    face_l2g = _partition_get_Yl2g_from_Xl2g(mesh.face2edge, edge_l2g)
    s_face2edge = _partition_get_local(mesh.face2edge, edge_g2l, face_l2g)

    face_g2l = _partition_index_l2g_to_g2l(face_l2g, nface(mesh))
    cell_l2g = _partition_get_Yl2g_from_Xl2g(mesh.cell2face, face_l2g)
    s_cell2face = _partition_get_local(mesh.cell2face, face_g2l, cell_l2g)

    return Mesh3D(s_vertices, s_edge2vert, s_face2edge, s_cell2face)
end

function submesh_e(mesh::Mesh3D, edge_l2g::AbstractVector{<:Integer})
    face_l2g = _partition_get_Yl2g_from_Xl2g(mesh.face2edge, edge_l2g)
    face_g2l = _partition_index_l2g_to_g2l(face_l2g, nface(mesh))
    cell_l2g = _partition_get_Yl2g_from_Xl2g(mesh.cell2face, face_l2g)
    s_cell2face = _partition_get_local(mesh.cell2face, face_g2l, cell_l2g)

    edge_g2l = _partition_index_l2g_to_g2l(edge_l2g, nedge(mesh))
    s_face2edge = _partition_get_local(mesh.face2edge, edge_g2l, face_l2g)

    vertex_l2g = _partition_get_Xl2g_from_Yl2g(mesh.edge2vert, edge_l2g, nvertex(mesh))
    vertex_g2l = _partition_index_l2g_to_g2l(vertex_l2g, nvertex(mesh))
    s_edge2vert = _partition_get_local(mesh.edge2vert, vertex_g2l, edge_l2g)

    s_vertices = mesh.vertex[:, vertex_l2g]
    return Mesh3D(s_vertices, s_edge2vert, s_face2edge, s_cell2face)
end

function submesh_f(mesh::Mesh3D, face_l2g::AbstractVector{<:Integer})
    cell_l2g = _partition_get_Yl2g_from_Xl2g(mesh.cell2face, face_l2g)
    face_g2l = _partition_index_l2g_to_g2l(face_l2g, nface(mesh))
    s_cell2face = _partition_get_local(mesh.cell2face, face_g2l, cell_l2g)

    edge_l2g = _partition_get_Xl2g_from_Yl2g(mesh.face2edge, face_l2g, nedge(mesh))
    edge_g2l = _partition_index_l2g_to_g2l(edge_l2g, nedge(mesh))
    s_face2edge = _partition_get_local(mesh.face2edge, edge_g2l, face_l2g)

    vertex_l2g = _partition_get_Xl2g_from_Yl2g(mesh.edge2vert, edge_l2g, nvertex(mesh))
    vertex_g2l = _partition_index_l2g_to_g2l(vertex_l2g, nvertex(mesh))
    s_edge2vert = _partition_get_local(mesh.edge2vert, vertex_g2l, edge_l2g)

    s_vertices = mesh.vertex[:, vertex_l2g]
    return Mesh3D(s_vertices, s_edge2vert, s_face2edge, s_cell2face)
end

function submesh_c(mesh::Mesh3D, cell_l2g::AbstractVector{<:Integer})
    face_l2g = _partition_get_Xl2g_from_Yl2g(mesh.cell2face, cell_l2g, nface(mesh))
    face_g2l = _partition_index_l2g_to_g2l(face_l2g, nface(mesh))
    s_cell2face = _partition_get_local(mesh.cell2face, face_g2l, cell_l2g)

    edge_l2g = _partition_get_Xl2g_from_Yl2g(mesh.face2edge, face_l2g, nedge(mesh))
    edge_g2l = _partition_index_l2g_to_g2l(edge_l2g, nedge(mesh))
    s_face2edge = _partition_get_local(mesh.face2edge, edge_g2l, face_l2g)

    vertex_l2g = _partition_get_Xl2g_from_Yl2g(mesh.edge2vert, edge_l2g, nvertex(mesh))
    vertex_g2l = _partition_index_l2g_to_g2l(vertex_l2g, nvertex(mesh))
    s_edge2vert = _partition_get_local(mesh.edge2vert, vertex_g2l, edge_l2g)

    s_vertices = mesh.vertex[:, vertex_l2g]
    return Mesh3D(s_vertices, s_edge2vert, s_face2edge, s_cell2face)
end

function submesh(mesh::Mesh3D, var::Symbol, l2g::AbstractVector{<:Integer})
    if var == :vertex
        return submesh_v(mesh, l2g)
    elseif var == :edge
        return submesh_e(mesh, l2g)
    elseif var == :face
        return submesh_f(mesh, l2g)
    elseif var == :cell
        return submesh_c(mesh, l2g)
    else
        error("var: $(var) not supported, availables are: :vertex, :edge, :face, :cell")
    end
end

function submesh(mesh::FEMesh; ks...)
    if length(ks) > 1
        error("too many keyword parameters")
    end
    k = keys(ks)[1]
    return submesh(mesh, k, ks[k])
end

#
#
#

export submeshdata

function submeshdata(mesh_l::Mesh1D, mesh_g::Mesh1D, data_g::Mesh1Ddata)
    vert_l2g = map(1:nvertex(mesh_l)) do ivl
        dists = map(eachindex(mesh_g.vertex)) do ivg
            abs(mesh_g.vertex[ivg] - mesh_l.vertex[ivl])
        end
        (vmin, imin) = findmin(dists)
        if vmin < 1e-5
            return imin
        else
            return 0
        end
    end
    if any(iszero, vert_l2g)
        error("any(iszero, vert_l2g)")
    end

    part_edge_g = _partition_get_global(mesh_l.edge2vert, vert_l2g)
    edge_l2g = map(axes(part_edge_g, 2)) do iel
        return _partition_match_col(_partition_match_vec_ignore_range,
            mesh_g.edge2vert, part_edge_g[:, iel])
    end

    dataidx = map(size(data_g.data), data_g.dimName) do n, k
        if k == :vertex
            return vert_l2g
        elseif k == :edge
            return edge_l2g
        else
            return n
        end
    end |> Tuple
    @info dataidx
    ndata = _partition_index_array(data_g.data, dataidx)

    return Mesh1Ddata(data_g.dimName, ndata)
end

function submeshdata(mesh_l::Mesh2D, mesh_g::Mesh2D, data_g::Mesh2Ddata)
    error("this function is being developed")
    return nothing
end

function submeshdata(mesh_l::Mesh3D, mesh_g::Mesh3D, data_g::Mesh3Ddata)
    error("this function is being developed")
    return nothing
end
