# See document Development.md for variable naming style
# simplify:
#   pI2tI -> I_p2t
#   tI2pI -> I_t2p
#   p_pIpJ2pK -> p_IJ2K
#   p_pJpIl2pI -> p_J2I

function _partition_It2p(I_p2t::Vector{Int}, t_ni::Int)
    I_t2p = zeros(Int, t_ni)
    foreach(eachindex(I_p2t)) do ipart
        itotal = I_p2t[ipart]
        I_t2p[itotal] = ipart
    end
    return I_t2p
end

function _partition_Ip2t(t_J2I::Matrix{Int}, J_p2t::Vector{Int})
    flag_tI = falses(maximum(t_J2I[:, J_p2t]))
    for tj = J_p2t
        for ti = t_J2I[:, tj]
            flag_tI[ti] = true
        end
    end
    return findall(flag_tI)
end

function _partition_Jp2t(t_J2I::Matrix{Int}, I_p2t::Vector{Int})
    flag_matched_tJ2I = map(ti -> ti in I_p2t, t_J2I)
    flag_tJ = all(flag_matched_tJ2I, dims=1) |> vec
    return findall(flag_tJ)
end

function _partition_p_J2I(t_J2I::Matrix{Int}, I_t2p::Vector{Int}, J_p2t::Vector{Int})
    return I_t2p[t_J2I[:, J_p2t]]
end

function _partition_t_pJ2tI(p_J2I::Matrix{Int}, I_p2t::Vector{Int})
    return I_p2t[p_J2I]
end

function _partition_colmatch(relation::Matrix{Int}, reference::Vector{Int})
    v0 = sort(reference)
    for ic in axes(relation, 2)
        v1 = sort(relation[:, ic])
        if v0 == v1
            return ic
        end
    end
    return 0
end

function _partition_Jp2t(t_J2I::Matrix{Int}, pJ2tI::Matrix{Int})
    return map(pj2tI -> _partition_colmatch(t_J2I, vec(pj2tI[:])), eachcol(pJ2tI))
end

const _APPROX_ZERO_THRESHOLD = 1e-5

function _partition_pv2tv(p_vert2coor::Vector{Float64}, t_vert2coor::Vector{Float64})
    return map(eachindex(p_vert2coor)) do pv
        dist = map(tv -> abs(p_vert2coor[pv] - t_vert2coor[tv]), eachindex(t_vert2coor))
        (vmin, imin) = findmin(dist)
        if vmin < _APPROX_ZERO_THRESHOLD
            return imin
        else
            return 0
        end
    end
end

function _partition_pv2tv(p_vert2coor::Matrix{Float64}, t_vert2coor::Matrix{Float64})
    return map(axes(p_vert2coor, 2)) do pv
        dist = map(tv -> norm(p_vert2coor[:, pv] - t_vert2coor[:, tv]), axes(t_vert2coor, 2))
        (vmin, imin) = findmin(dist)
        if vmin < _APPROX_ZERO_THRESHOLD
            return imin
        else
            return 0
        end
    end
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

#
#
#

export submesh

"""
```julia
submesh(tmesh::FEMesh, var::Symbol, part2totl::AbstractVector{<:Integer}) -> FEMesh
```
"""
function submesh end

function submesh_v(tmesh::Mesh1D, vert_p2t::AbstractVector{<:Integer})
    vert_t2p = _partition_It2p(vert_p2t, nvertex(tmesh))
    edge_p2t = _partition_Jp2t(tmesh.edge2vert, vert_p2t)

    p_vert2coor = tmesh.vert2coor[:, vert_p2t]
    p_edge2vert = _partition_p_J2I(tmesh.edge2vert, vert_t2p, edge_p2t)
    return Mesh1D(p_vert2coor, p_edge2vert)
end

function submesh_e(tmesh::Mesh1D, edge_p2t::AbstractVector{<:Integer})
    vert_p2t = _partition_Ip2t(tmesh.edge2vert, edge_p2t)
    vert_t2p = _partition_It2p(vert_p2t, nvertex(tmesh))

    p_vert2coor = tmesh.vert2coor[:, vert_p2t]
    p_edge2vert = _partition_p_J2I(tmesh.edge2vert, vert_t2p, edge_p2t)
    return Mesh1D(p_vert2coor, p_edge2vert)
end

function submesh(tmesh::Mesh1D, var::Symbol, part2totl::AbstractVector{<:Integer})
    if var == :vertex
        return submesh_v(tmesh, part2totl)
    elseif var == :edge
        return submesh_e(tmesh, part2totl)
    else
        error("var: $(var) not supported, availables are: :vertex, :edge")
    end
end

function submesh_v(tmesh::Mesh2D, vert_p2t::AbstractVector{<:Integer})
    vert_t2p = _partition_It2p(vert_p2t, nvertex(tmesh))
    edge_p2t = _partition_Jp2t(tmesh.edge2vert, vert_p2t)
    edge_t2p = _partition_It2p(edge_p2t, nedge(tmesh))
    face_p2t = _partition_Jp2t(tmesh.face2edge, edge_p2t)

    p_vert2coor = tmesh.vert2coor[:, vert_p2t]
    p_edge2vert = _partition_p_J2I(tmesh.edge2vert, vert_t2p, edge_p2t)
    p_face2edge = _partition_p_J2I(tmesh.face2edge, edge_t2p, face_p2t)
    return Mesh2D(p_vert2coor, p_edge2vert, p_face2edge)
end

function submesh_e(tmesh::Mesh2D, edge_p2t::AbstractVector{<:Integer})
    vert_p2t = _partition_Ip2t(tmesh.edge2vert, edge_p2t)
    vert_t2p = _partition_It2p(vert_p2t, nvertex(tmesh))
    edge_t2p = _partition_It2p(edge_p2t, nedge(tmesh))
    face_p2t = _partition_Jp2t(tmesh.face2edge, edge_p2t)

    p_vert2coor = tmesh.vert2coor[:, vert_p2t]
    p_edge2vert = _partition_p_J2I(tmesh.edge2vert, vert_t2p, edge_p2t)
    p_face2edge = _partition_p_J2I(tmesh.face2edge, edge_t2p, face_p2t)
    return Mesh2D(p_vert2coor, p_edge2vert, p_face2edge)
end

function submesh_f(tmesh::Mesh2D, face_p2t::AbstractVector{<:Integer})
    edge_p2t = _partition_Ip2t(tmesh.face2edge, face_p2t)
    edge_t2p = _partition_It2p(edge_p2t, nedge(tmesh))
    vert_p2t = _partition_Ip2t(tmesh.edge2vert, edge_p2t)
    vert_t2p = _partition_It2p(vert_p2t, nvertex(tmesh))

    p_vert2coor = tmesh.vert2coor[:, vert_p2t]
    p_edge2vert = _partition_p_J2I(tmesh.edge2vert, vert_t2p, edge_p2t)
    p_face2edge = _partition_p_J2I(tmesh.face2edge, edge_t2p, face_p2t)
    return Mesh2D(p_vert2coor, p_edge2vert, p_face2edge)
end

function submesh(tmesh::Mesh2D, var::Symbol, part2totl::AbstractVector{<:Integer})
    if var == :vertex
        return submesh_v(tmesh, part2totl)
    elseif var == :edge
        return submesh_e(tmesh, part2totl)
    elseif var == :face
        return submesh_f(tmesh, part2totl)
    else
        error("var: $(var) not supported, availables are: :vertex, :edge, :face")
    end
end

function submesh_v(tmesh::Mesh3D, vert_p2t::AbstractVector{<:Integer})
    vert_t2p = _partition_It2p(vert_p2t, nvertex(tmesh))
    edge_p2t = _partition_Jp2t(tmesh.edge2vert, vert_p2t)
    edge_t2p = _partition_It2p(edge_p2t, nedge(tmesh))
    face_p2t = _partition_Jp2t(tmesh.face2edge, edge_p2t)
    face_t2p = _partition_It2p(face_p2t, nface(tmesh))
    cell_p2t = _partition_Jp2t(tmesh.cell2face, face_p2t)

    p_vert2coor = tmesh.vert2coor[:, vert_p2t]
    p_edge2vert = _partition_p_J2I(tmesh.edge2vert, vert_t2p, edge_p2t)
    p_face2edge = _partition_p_J2I(tmesh.face2edge, edge_t2p, face_p2t)
    p_cell2face = _partition_p_J2I(tmesh.cell2face, face_t2p, cell_p2t)
    return Mesh3D(p_vert2coor, p_edge2vert, p_face2edge, p_cell2face)
end

function submesh_e(tmesh::Mesh3D, edge_p2t::AbstractVector{<:Integer})
    vert_p2t = _partition_Ip2t(mesh.edge2vert, edge_p2t)
    vert_t2p = _partition_It2p(vert_p2t, nvertex(tmesh))
    edge_t2p = _partition_It2p(edge_p2t, nedge(tmesh))
    face_p2t = _partition_Jp2t(tmesh.face2edge, edge_p2t)
    face_t2p = _partition_It2p(face_p2t, nface(tmesh))
    cell_p2t = _partition_Jp2t(tmesh.cell2face, face_p2t)

    p_vert2coor = tmesh.vert2coor[:, vert_p2t]
    p_edge2vert = _partition_p_J2I(tmesh.edge2vert, vert_t2p, edge_p2t)
    p_face2edge = _partition_p_J2I(tmesh.face2edge, edge_t2p, face_p2t)
    p_cell2face = _partition_p_J2I(tmesh.cell2face, face_t2p, cell_p2t)
    return Mesh3D(p_vert2coor, p_edge2vert, p_face2edge, p_cell2face)
end

function submesh_f(tmesh::Mesh3D, face_p2t::AbstractVector{<:Integer})
    edge_p2t = _partition_Ip2t(tmesh.face2edge, face_p2t)
    edge_t2p = _partition_It2p(edge_p2t, nedge(tmesh))
    vert_p2t = _partition_Ip2t(tmesh.edge2vert, edge_p2t)
    vert_t2p = _partition_It2p(vert_p2t, nvertex(tmesh))
    face_t2p = _partition_It2p(face_p2t, nface(tmesh))
    cell_p2t = _partition_Jp2t(tmesh.cell2face, face_p2t)

    p_vert2coor = tmesh.vert2coor[:, vert_p2t]
    p_edge2vert = _partition_p_J2I(tmesh.edge2vert, vert_t2p, edge_p2t)
    p_face2edge = _partition_p_J2I(tmesh.face2edge, edge_t2p, face_p2t)
    p_cell2face = _partition_p_J2I(tmesh.cell2face, face_t2p, cell_p2t)
    return Mesh3D(p_vert2coor, p_edge2vert, p_face2edge, p_cell2face)
end

function submesh_c(tmesh::Mesh3D, cell_p2t::AbstractVector{<:Integer})
    face_p2t = _partition_Ip2t(tmesh.cell2face, cell_p2t)
    face_t2p = _partition_It2p(face_p2t, nface(tmesh))
    edge_p2t = _partition_Ip2t(tmesh.face2edge, face_p2t)
    edge_t2p = _partition_It2p(edge_p2t, nedge(tmesh))
    vert_p2t = _partition_Ip2t(tmesh.edge2vert, edge_p2t)
    vert_t2p = _partition_It2p(vert_p2t, nvertex(tmesh))

    p_vert2coor = tmesh.vert2coor[:, vert_p2t]
    p_edge2vert = _partition_p_J2I(tmesh.edge2vert, vert_t2p, edge_p2t)
    p_face2edge = _partition_p_J2I(tmesh.face2edge, edge_t2p, face_p2t)
    p_cell2face = _partition_p_J2I(tmesh.cell2face, face_t2p, cell_p2t)
    return Mesh3D(p_vert2coor, p_edge2vert, p_face2edge, p_cell2face)
end

function submesh(tmesh::Mesh3D, var::Symbol, part2totl::AbstractVector{<:Integer})
    if var == :vertex
        return submesh_v(tmesh, part2totl)
    elseif var == :edge
        return submesh_e(tmesh, part2totl)
    elseif var == :face
        return submesh_f(tmesh, part2totl)
    elseif var == :cell
        return submesh_c(tmesh, part2totl)
    else
        error("var: $(var) not supported, availables are: :vertex, :edge, :face, :cell")
    end
end

"""
```julia
submesh(tmesh::FEMesh; vertex/edge/face/cell) -> FEMesh
```
"""
function submesh(tmesh::FEMesh; ks...)
    if length(ks) > 1
        error("too many keyword parameters")
    end
    k = keys(ks)[1]
    return submesh(tmesh, k, ks[k])
end

#
#
#

export submeshdata

"""
```julia
submeshdata(partmesh::FEMesh, totalmesh::FEMesh, t_data::Meshdata) -> Meshdata
```
"""
function submeshdata end

function submeshdata(pmesh::Mesh1D, tmesh::Mesh1D, t_data::Mesh1Ddata)
    vert_p2t = _partition_pv2tv(pmesh.vert2coor, tmesh.vert2coor)
    if any(iszero, vert_p2t) || isempty(vert_p2t)
        error("any(iszero, vert_l2g)")
    end

    pedge2tvert = _partition_t_pJ2tI(pmesh.edge2vert, vert_p2t)
    edge_p2t = _partition_Jp2t(tmesh.edge2vert, pedge2tvert)

    dataidx = map(size(t_data.data), t_data.dimName) do n, k
        if k == :vertex
            return vert_p2t
        elseif k == :edge
            return edge_p2t
        else
            return n
        end
    end |> Tuple
    @info dataidx
    ndata = _partition_index_array(t_data.data, dataidx)

    return Mesh1Ddata(t_data.dimName, ndata)
end

function submeshdata(pmesh::Mesh2D, tmesh::Mesh2D, t_data::Mesh2Ddata)
    vert_p2t = _partition_pv2tv(pmesh.vertex, tmesh.vertex)
    if any(iszero, vert_l2g) || isempty(vert_p2t)
        error("any(iszero, vert_l2g)")
    end
    pedge2tvert = _partition_t_pJ2tI(pmesh.edge2vert, vert_p2t)
    edge_p2t = _partition_Jp2t(tmesh.edge2vert, pedge2tvert)
    pface2tedge = _partition_t_pJ2tI(pmesh.face2edge, edge_p2t)
    face_p2t = _partition_Jp2t(tmesh.face2edge, pface2tedge)

    dataidx = map(size(t_data.data), t_data.dimName) do n, k
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
    ndata = _partition_index_array(t_data.data, dataidx)

    return Mesh2Ddata(t_data.dimName, ndata)
end

function submeshdata(pmesh::Mesh3D, tmesh::Mesh3D, t_data::Mesh3Ddata)
    vert_p2t = _partition_pv2tv(pmesh.vertex, tmesh.vertex)
    if any(iszero, vert_l2g) || isempty(vert_p2t)
        error("any(iszero, vert_l2g)")
    end
    pedge2tvert = _partition_t_pJ2tI(pmesh.edge2vert, vert_p2t)
    edge_p2t = _partition_Jp2t(tmesh.edge2vert, pedge2tvert)
    pface2tedge = _partition_t_pJ2tI(pmesh.face2edge, edge_p2t)
    face_p2t = _partition_Jp2t(tmesh.face2edge, pface2tedge)
    pcell2tface = _partition_t_pJ2tI(pmesh.cell2face, face_p2t)
    cell_p2t = _partition_Jp2t(tmesh.cell2face, pcell2tface)

    dataidx = map(size(t_data.data), t_data.dimName) do n, k
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
    ndata = _partition_index_array(t_data.data, dataidx)

    return Mesh3Ddata(t_data.dimName, ndata)
end
