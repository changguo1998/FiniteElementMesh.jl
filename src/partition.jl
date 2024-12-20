# To make differences of same data between part(new) mesh and total(old) mesh,
# the variable is appended with suffix `_p` and `_t`
# suffix `_p2t` and `_t2p` is also added to show the usage for index converting
# suffix `_l`, `_g`, `_l2g` and `_g2l` is used to distinguish index in *same* mesh

# TODO check and replace _l, _g suffix
# TODO current progress:

function _partition_index_p2t_to_t2p(index_p2t::Vector{Int}, n_total::Int)
    index_t2p = zeros(Int, n_total)
    foreach(eachindex(index_p2t)) do ipart
        itotal = index_p2t[ipart]
        index_t2p[itotal] = ipart
    end
    return index_t2p
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

const _APPROX_ZERO_THRESHOLD = 1e-5

_partition_match_vec(vec1::Vector{Int}, vec2::Vector{Int}) = vec1 == vec2
_partition_match_vec(vec1::Vector{Float64}, vec2::Vector{Float64}) = vec1 == vec2

function _partition_match_vec_ignore_range(vec1::Vector{Int}, vec2::Vector{Int})
    v1 = sort(vec1)
    v2 = sort(vec2)
    return v1 == v2
end

function _partition_match_col_idx(relation::Matrix{Int}, reference::Vector{Int})
    for ic in axes(relation, 2)
        if _partition_match_vec_ignore_range(relation[:, ic], reference)
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

function _partition_vertex_l2g_from_coordinate(vertex_p::Vector{Float64}, vert_t::Vector{Float64})
    return map(eachindex(vertex_p)) do ivert_p
        dist = map(ivert_t->abs(vertex_p[ivert_p]-vert_t[ivert_t]), eachindex(vert_t))
        (vmin, imin) = findmin(dist)
        if vmin < _APPROX_ZERO_THRESHOLD
            return imin
        else
            return 0
        end
    end
end

function _partition_vertex_l2g_from_coordinate(vertex_p::Matrix{Float64}, vert_t::Matrix{Float64})
    return map(axes(vertex_p, 2)) do ivert_p
        dist = map(ivert_t->norm(vertex_p[:, ivert_p]-vert_t[:, ivert_t]), axes(vert_t, 2))
        (vmin, imin) = findmin(dist)
        if vmin < _APPROX_ZERO_THRESHOLD
            return imin
        else
            return 0
        end
    end
end

function _partition_get_l2g_from_index(part_Y2X_g::Matrix{Int}, Y2X_g::Matrix{Int})
    return map(col_part_l->_partition_match_col_idx(Y2X_g, col_part_l), eachcol(part_Y2X_g))
end

#
#
#

export submesh

function submesh_v(mesh::Mesh1D, vertex_l2g::AbstractVector{<:Integer})
    s_vertices = mesh.vertex[vertex_l2g]

    vertex_g2l = _partition_index_p2t_to_t2p(vertex_l2g, nvertex(mesh))
    edge_l2g = _partition_get_Yl2g_from_Xl2g(mesh.edge2vert, vertex_l2g)
    s_edge2vert = _partition_get_local(mesh.edge2vert, vertex_g2l, edge_l2g)
    return Mesh1D(s_vertices, s_edge2vert)
end

function submesh_e(mesh::Mesh1D, edge_l2g::AbstractVector{<:Integer})
    vertex_l2g = _partition_get_Xl2g_from_Yl2g(mesh.edge2vert, edge_l2g, nvertex(mesh))
    vertex_g2l = _partition_index_p2t_to_t2p(vertex_l2g, nvertex(mesh))

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

    vertex_g2l = _partition_index_p2t_to_t2p(vertex_l2g, nvertex(mesh))
    edge_l2g = _partition_get_Yl2g_from_Xl2g(mesh.edge2vert, vertex_l2g)
    s_edge2vert = _partition_get_local(mesh.edge2vert, vertex_g2l, edge_l2g)

    edge_g2l = _partition_index_p2t_to_t2p(edge_l2g, nedge(mesh))
    face_l2g = _partition_get_Yl2g_from_Xl2g(mesh.face2edge, edge_l2g)
    s_face2edge = _partition_get_local(mesh.face2edge, edge_g2l, face_l2g)
    return Mesh2D(s_vertices, s_edge2vert, s_face2edge)
end

function submesh_e(mesh::Mesh2D, edge_l2g::AbstractVector{<:Integer})
    face_l2g = _partition_get_Yl2g_from_Xl2g(mesh.face2edge, edge_l2g)
    edge_g2l = _partition_index_p2t_to_t2p(edge_l2g, nedge(mesh))
    s_face2edge = _partition_get_local(mesh.face2edge, edge_g2l, face_l2g)

    vertex_l2g = _partition_get_Xl2g_from_Yl2g(mesh.edge2vert, edge_l2g, nvertex(mesh))
    vertex_g2l = _partition_index_p2t_to_t2p(vertex_l2g, nvertex(mesh))
    s_edge2vert = _partition_get_local(mesh.edge2vert, vertex_g2l, edge_l2g)

    s_vertices = mesh.vertex[:, vertex_l2g]
    return Mesh2D(s_vertices, s_edge2vert, s_face2edge)
end

function submesh_f(mesh::Mesh2D, face_l2g::AbstractVector{<:Integer})
    edge_l2g = _partition_get_Xl2g_from_Yl2g(mesh.face2edge, face_l2g, nedge(mesh))
    edge_g2l = _partition_index_p2t_to_t2p(edge_l2g, nedge(mesh))
    s_face2edge = _partition_get_local(mesh.face2edge, edge_g2l, face_l2g)

    vertex_l2g = _partition_get_Xl2g_from_Yl2g(mesh.edge2vert, edge_l2g, nvertex(mesh))
    vertex_g2l = _partition_index_p2t_to_t2p(vertex_l2g, nvertex(mesh))
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

    vertex_g2l = _partition_index_p2t_to_t2p(vertex_l2g, nvertex(mesh))
    edge_l2g = _partition_get_Yl2g_from_Xl2g(mesh.edge2vert, vertex_l2g)
    s_edge2vert = _partition_get_local(mesh.edge2vert, vertex_g2l, edge_l2g)

    edge_g2l = _partition_index_p2t_to_t2p(edge_l2g, nedge(mesh))
    face_l2g = _partition_get_Yl2g_from_Xl2g(mesh.face2edge, edge_l2g)
    s_face2edge = _partition_get_local(mesh.face2edge, edge_g2l, face_l2g)

    face_g2l = _partition_index_p2t_to_t2p(face_l2g, nface(mesh))
    cell_l2g = _partition_get_Yl2g_from_Xl2g(mesh.cell2face, face_l2g)
    s_cell2face = _partition_get_local(mesh.cell2face, face_g2l, cell_l2g)

    return Mesh3D(s_vertices, s_edge2vert, s_face2edge, s_cell2face)
end

function submesh_e(mesh::Mesh3D, edge_l2g::AbstractVector{<:Integer})
    face_l2g = _partition_get_Yl2g_from_Xl2g(mesh.face2edge, edge_l2g)
    face_g2l = _partition_index_p2t_to_t2p(face_l2g, nface(mesh))
    cell_l2g = _partition_get_Yl2g_from_Xl2g(mesh.cell2face, face_l2g)
    s_cell2face = _partition_get_local(mesh.cell2face, face_g2l, cell_l2g)

    edge_g2l = _partition_index_p2t_to_t2p(edge_l2g, nedge(mesh))
    s_face2edge = _partition_get_local(mesh.face2edge, edge_g2l, face_l2g)

    vertex_l2g = _partition_get_Xl2g_from_Yl2g(mesh.edge2vert, edge_l2g, nvertex(mesh))
    vertex_g2l = _partition_index_p2t_to_t2p(vertex_l2g, nvertex(mesh))
    s_edge2vert = _partition_get_local(mesh.edge2vert, vertex_g2l, edge_l2g)

    s_vertices = mesh.vertex[:, vertex_l2g]
    return Mesh3D(s_vertices, s_edge2vert, s_face2edge, s_cell2face)
end

function submesh_f(mesh::Mesh3D, face_l2g::AbstractVector{<:Integer})
    cell_l2g = _partition_get_Yl2g_from_Xl2g(mesh.cell2face, face_l2g)
    face_g2l = _partition_index_p2t_to_t2p(face_l2g, nface(mesh))
    s_cell2face = _partition_get_local(mesh.cell2face, face_g2l, cell_l2g)

    edge_l2g = _partition_get_Xl2g_from_Yl2g(mesh.face2edge, face_l2g, nedge(mesh))
    edge_g2l = _partition_index_p2t_to_t2p(edge_l2g, nedge(mesh))
    s_face2edge = _partition_get_local(mesh.face2edge, edge_g2l, face_l2g)

    vertex_l2g = _partition_get_Xl2g_from_Yl2g(mesh.edge2vert, edge_l2g, nvertex(mesh))
    vertex_g2l = _partition_index_p2t_to_t2p(vertex_l2g, nvertex(mesh))
    s_edge2vert = _partition_get_local(mesh.edge2vert, vertex_g2l, edge_l2g)

    s_vertices = mesh.vertex[:, vertex_l2g]
    return Mesh3D(s_vertices, s_edge2vert, s_face2edge, s_cell2face)
end

function submesh_c(mesh::Mesh3D, cell_l2g::AbstractVector{<:Integer})
    face_l2g = _partition_get_Xl2g_from_Yl2g(mesh.cell2face, cell_l2g, nface(mesh))
    face_g2l = _partition_index_p2t_to_t2p(face_l2g, nface(mesh))
    s_cell2face = _partition_get_local(mesh.cell2face, face_g2l, cell_l2g)

    edge_l2g = _partition_get_Xl2g_from_Yl2g(mesh.face2edge, face_l2g, nedge(mesh))
    edge_g2l = _partition_index_p2t_to_t2p(edge_l2g, nedge(mesh))
    s_face2edge = _partition_get_local(mesh.face2edge, edge_g2l, face_l2g)

    vertex_l2g = _partition_get_Xl2g_from_Yl2g(mesh.edge2vert, edge_l2g, nvertex(mesh))
    vertex_g2l = _partition_index_p2t_to_t2p(vertex_l2g, nvertex(mesh))
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

function submeshdata(mesh_part::Mesh1D, mesh_total::Mesh1D, data_total::Mesh1Ddata)
    vert_p2t = _partition_vertex_l2g_from_coordinate(mesh_part.vertex, mesh_total.vertex)
    if any(iszero, vert_p2t) || isempty(vert_p2t)
        error("any(iszero, vert_l2g)")
    end

    part_edge2vert_g = _partition_get_global(mesh_part.edge2vert, vert_p2t)
    edge_p2t = _partition_get_l2g_from_index(part_edge2vert_g, mesh_total.edge2vert)

    dataidx = map(size(data_total.data), data_total.dimName) do n, k
        if k == :vertex
            return vert_p2t
        elseif k == :edge
            return edge_p2t
        else
            return n
        end
    end |> Tuple
    @info dataidx
    ndata = _partition_index_array(data_total.data, dataidx)

    return Mesh1Ddata(data_total.dimName, ndata)
end

function submeshdata(mesh_part::Mesh2D, mesh_total::Mesh2D, data_total::Mesh2Ddata)
    vert_p2t = _partition_vertex_l2g_from_coordinate(mesh_part.vertex, mesh_total.vertex)
    if any(iszero, vert_l2g)
        error("any(iszero, vert_l2g)")
    end

    part_edge2vert_g = _partition_get_global(mesh_part.edge2vert, vert_p2t)
    edge_p2t = _partition_get_l2g_from_index(part_edge2vert_g, mesh_total.edge2vert)

    part_face2edge_g = _partition_get_global(mesh_part.face2edge, edge_p2t)
    face_p2t = _partition_get_l2g_from_index(part_face2edge_g, mesh_total.face2edge)


    dataidx = map(size(data_total.data), data_total.dimName) do n, k
        if k == :vertex
            return vert_p2t
        elseif k == :edge
            return edge_p2t
        elseif k == :face
            return face_p2t
        else
            return n
        end
    end |> Tuple
    @info dataidx
    ndata = _partition_index_array(data_total.data, dataidx)

    return Mesh2Ddata(data_total.dimName, ndata)
end

function submeshdata(mesh_part::Mesh3D, mesh_total::Mesh3D, data_total::Mesh3Ddata)
    error("this function is being developed")
    vert_p2t = _partition_vertex_l2g_from_coordinate(mesh_part.vertex, mesh_total.vertex)
    if any(iszero, vert_l2g)
        error("any(iszero, vert_l2g)")
    end

    part_edge2vert_g = _partition_get_global(mesh_part.edge2vert, vert_p2t)
    edge_p2t = _partition_get_l2g_from_index(part_edge2vert_g, mesh_total.edge2vert)

    part_face2edge_g = _partition_get_global(mesh_part.face2edge, edge_p2t)
    face_p2t = _partition_get_l2g_from_index(part_face2edge_g, mesh_total.face2edge)

    part_cell2face_g = _partition_get_global(mesh_part.cell2face, face_p2t)
    cell_p2t = _partition_get_l2g_from_index(part_cell2face_g, mesh_total.cell2face)

    dataidx = map(size(data_total.data), data_total.dimName) do n, k
        if k == :vertex
            return vert_p2t
        elseif k == :edge
            return edge_p2t
        elseif k == :face
            return face_p2t
        elseif k == :cell
            return cell_p2t
        else
            return n
        end
    end |> Tuple
    @info dataidx
    ndata = _partition_index_array(data_total.data, dataidx)

    return Mesh3Ddata(data_total.dimName, ndata)
end
