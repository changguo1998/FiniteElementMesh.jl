
#
# 0. basic function
#

function _relation_high2low_skip_level(Z2Y::Matrix{Int}, Y2X::Matrix{Int})
    nz = size(Z2Y, 2)
    # ny = size(Y2X, 2)
    ny_per_z = size(Z2Y, 1)
    nx_per_y = size(Y2X, 1)

    buffer = zeros(Int, ny_per_z * nx_per_y, nz)
    for iz = axes(buffer, 2)
        for iyl = axes(Z2Y, 1)
            yshift = (iyl - 1) * nx_per_y
            for ix = axes(Y2X, 1)
                buffer[ix+yshift, iz] = Y2X[ix, Z2Y[iyl, iz]]
            end
        end
    end
    nx_per_z = maximum(map(col -> length(unique(col)), eachcol(buffer)))
    Z2X = zeros(Int, nx_per_z, nz)
    for iz = axes(Z2X, 2)
        Z2X[:, iz] .= unique(buffer[:, iz])
    end
    return Z2X
end

function _relation_get_low2high_same_level(Y2X::Matrix{Int}, NX::Integer)
    ncountx = zeros(Int, NX)
    for ix = Y2X
        ncountx[ix] += 1
    end
    ny_per_x = maximum(ncountx)
    X2Y = zeros(Int, ny_per_x, NX)
    for iy = axes(Y2X, 2)
        for ixl = axes(Y2X, 1)
            ix = Y2X[ixl, iy]
            if iy in X2Y[:, ix]
                continue
            end
            for iyl = axes(X2Y, 1)
                if iszero(X2Y[iyl, ix])
                    X2Y[iyl, ix] = iy
                    break
                end
            end
        end
    end
    return X2Y
end

function _relation_get_low2high_id_same_level(Y2X::Matrix{Int}, X2Y::Matrix{Int})
    X2Yl = zeros(Int, size(X2Y))
    for ix = axes(X2Y, 2), iyl = axes(X2Y, 1)
        iy = X2Y[iyl, ix]
        if iszero(iy)
            continue
        end
        for ixl = axes(Y2X, 1)
            if Y2X[ixl, iy] == ix
                X2Yl[iyl, ix] = ixl
            end
        end
    end
    return X2Yl
end

function _relation_low2high_skip_level(X2Y::Matrix{Int}, Y2Z::Matrix{Int})
    nx = size(X2Y, 2)
    ny_per_x = size(X2Y, 1)
    nz_per_y = size(Y2Z, 1)

    buffer = zeros(Int, ny_per_x * nz_per_y, nx)
    for ix = axes(buffer, 2)
        for iyl = axes(X2Y, 1)
            if iszero(X2Y[iyl, ix])
                continue
            end
            yshift = (iyl - 1) * nz_per_y
            for izl = axes(Y2Z, 1)
                buffer[izl+yshift, ix] = Y2Z[izl, X2Y[iyl, ix]]
            end
        end
    end
    nz_per_x = maximum(map(col -> length(unique(filter(!iszero, col))), eachcol(buffer)))
    X2Z = zeros(Int, nz_per_x, nx)
    for ix = axes(X2Z, 2)
        t = unique(filter(!iszero, buffer[:, ix]))
        X2Z[eachindex(t), ix] .= t
    end
    return X2Z
end

#
# 1. auxillary_relation
#

export auxillary_relation

"""
```
auxillary_relation(mesh::Mesh1D) -> Mesh1DAux
```
"""
function auxillary_relation(mesh::Mesh1D)
    v2e = _relation_get_low2high_same_level(mesh.edge2vert, nvertex(mesh))
    v2e_l = _relation_get_low2high_id_same_level(mesh.edge2vert, v2e)
    return Mesh1Daux(v2e, v2e_l)
end

"""
```
auxillary_relation(mesh::Mesh2D) -> Mesh2DAux
```
"""
function auxillary_relation(mesh::Mesh2D)
    f2v = _relation_high2low_skip_level(mesh.face2edge, mesh.edge2vert)

    e2f = _relation_get_low2high_same_level(mesh.face2edge, nedge(mesh))
    v2e = _relation_get_low2high_same_level(mesh.edge2vert, nvertex(mesh))
    v2f = _relation_low2high_skip_level(v2e, e2f)

    e2f_l = _relation_get_low2high_id_same_level(mesh.face2edge, e2f)
    v2e_l = _relation_get_low2high_id_same_level(mesh.edge2vert, v2e)
    v2f_l = _relation_get_low2high_id_same_level(f2v, v2f)
    return Mesh2Daux(f2v, v2e, v2f, e2f, v2e_l, v2f_l, e2f_l)
end

"""
```
auxillary_relation(mesh::Mesh3D) -> Mesh3DAux
```
"""
function auxillary_relation(mesh::Mesh3D)
    f2v = _relation_high2low_skip_level(mesh.face2edge, mesh.edge2vert)
    c2e = _relation_high2low_skip_level(mesh.cell2face, mesh.face2edge)
    c2v = _relation_high2low_skip_level(mesh.cell2face, f2v)

    v2e = _relation_get_low2high_same_level(mesh.edge2vert, nvertex(mesh))
    e2f = _relation_get_low2high_same_level(mesh.face2edge, nedge(mesh))
    f2c = _relation_get_low2high_same_level(mesh.cell2face, nface(mesh))
    v2f = _relation_low2high_skip_level(v2e, e2f)
    e2c = _relation_low2high_skip_level(e2f, f2c)
    v2c = _relation_low2high_skip_level(v2f, f2c)

    v2e_l = _relation_get_low2high_id_same_level(mesh.edge2vert, v2e)
    e2f_l = _relation_get_low2high_id_same_level(mesh.face2edge, e2f)
    f2c_l = _relation_get_low2high_id_same_level(mesh.cell2face, f2c)
    v2f_l = _relation_get_low2high_id_same_level(f2v, v2f)
    e2c_l = _relation_get_low2high_id_same_level(c2e, e2c)
    v2c_l = _relation_get_low2high_id_same_level(c2v, v2c)
    return Mesh3Daux(f2v, c2e, c2v, v2e, v2f, v2c, e2f, e2c, f2c, v2e_l, v2f_l, v2c_l, e2f_l, e2c_l, f2c_l)
end
