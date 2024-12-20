
#
# 0. basic function
#

function _relation_high2low_transfer(K2J::Matrix{Int}, J2I::Matrix{Int})
    nk = size(K2J, 2)
    # ny = size(Y2X, 2)
    nj_per_k = size(K2J, 1)
    ni_per_j = size(J2I, 1)

    buffer = zeros(Int, nj_per_k * ni_per_j, nk)
    ils = axes(J2I, 1)
    for k = axes(buffer, 2), jl = axes(K2J, 1)
        j = K2J[jl, k]
        jshift = (jl - 1) * ni_per_j
        buffer[ils.+jshift, k] .= J2I[:, j]
        # for il = axes(J2I, 1)
        #     buffer[il+jshift, k] = J2I[il, j]
        # end
    end
    ni_per_k = maximum(map(col -> length(unique(col)), eachcol(buffer)))
    K2I = zeros(Int, ni_per_k, nk)
    for k = axes(K2I, 2)
        K2I[:, k] .= unique(buffer[:, k])
    end
    return K2I
end

function _relation_low2high_transfer(I2J::Matrix{Int}, J2K::Matrix{Int})
    ni = size(I2J, 2)
    nj_per_i = size(I2J, 1)
    nk_per_j = size(J2K, 1)

    buffer = zeros(Int, nj_per_i * nk_per_j, ni)
    kls = axes(J2K, 1)
    for i = axes(buffer, 2), jl = axes(I2J, 1)
        j = I2J[jl, i]
        if iszero(j)
            continue
        end
        jshift = (jl - 1) * nk_per_j
        buffer[kls .+ jshift, i] = J2K[:, j]
        # for kl = axes(J2K, 1)
        #     buffer[kl+jshift, i] = J2K[kl, j]
        # end
    end
    nk_per_i = maximum(map(col -> length(unique(filter(!iszero, col))), eachcol(buffer)))
    I2K = zeros(Int, nk_per_i, ni)
    for i = axes(I2K, 2)
        t = unique(filter(!iszero, buffer[:, i]))
        I2K[eachindex(t), i] .= t
    end
    return I2K
end

function _relation_I2J(J2I::Matrix{Int}, ni::Integer)
    ncountx = zeros(Int, ni)
    for i = J2I
        ncountx[i] += 1
    end
    nj_per_i = maximum(ncountx)
    I2J = zeros(Int, nj_per_i, ni)
    for j = axes(J2I, 2), il = axes(J2I, 1)
        i = J2I[il, j]
        if j in I2J[:, i]
            continue
        end
        for jl = axes(I2J, 1)
            if iszero(I2J[jl, i])
                I2J[jl, i] = j
                break
            end
        end
    end
    return I2J
end

function _relation_IgJl2Il(J2I::Matrix{Int}, I2J::Matrix{Int})
    IgJl2Il = zeros(Int, size(I2J))
    for i = axes(I2J, 2), jl = axes(I2J, 1)
        j = I2J[jl, i]
        if iszero(j)
            continue
        end
        for il = axes(J2I, 1)
            if J2I[il, j] == i
                IgJl2Il[jl, i] = il
            end
        end
    end
    return IgJl2Il
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
    v2e = _relation_I2J(mesh.edge2vert, nvertex(mesh))
    vel2vl = _relation_IgJl2Il(mesh.edge2vert, v2e)
    return Mesh1Daux(v2e, vel2vl)
end

"""
```
auxillary_relation(mesh::Mesh2D) -> Mesh2DAux
```
"""
function auxillary_relation(mesh::Mesh2D)
    f2v = _relation_high2low_transfer(mesh.face2edge, mesh.edge2vert)

    e2f = _relation_I2J(mesh.face2edge, nedge(mesh))
    v2e = _relation_I2J(mesh.edge2vert, nvertex(mesh))
    v2f = _relation_low2high_transfer(v2e, e2f)

    efl2el = _relation_IgJl2Il(mesh.face2edge, e2f)
    vel2vl = _relation_IgJl2Il(mesh.edge2vert, v2e)
    vfl2vl = _relation_IgJl2Il(f2v, v2f)
    return Mesh2Daux(f2v, v2e, v2f, e2f, vel2vl, vfl2vl, efl2el)
end

"""
```
auxillary_relation(mesh::Mesh3D) -> Mesh3DAux
```
"""
function auxillary_relation(mesh::Mesh3D)
    f2v = _relation_high2low_transfer(mesh.face2edge, mesh.edge2vert)
    c2e = _relation_high2low_transfer(mesh.cell2face, mesh.face2edge)
    c2v = _relation_high2low_transfer(mesh.cell2face, f2v)

    v2e = _relation_I2J(mesh.edge2vert, nvertex(mesh))
    e2f = _relation_I2J(mesh.face2edge, nedge(mesh))
    f2c = _relation_I2J(mesh.cell2face, nface(mesh))
    v2f = _relation_low2high_transfer(v2e, e2f)
    e2c = _relation_low2high_transfer(e2f, f2c)
    v2c = _relation_low2high_transfer(v2f, f2c)

    vel2vl = _relation_IgJl2Il(mesh.edge2vert, v2e)
    efl2el = _relation_IgJl2Il(mesh.face2edge, e2f)
    fcl2fl = _relation_IgJl2Il(mesh.cell2face, f2c)
    vfl2vl = _relation_IgJl2Il(f2v, v2f)
    ecl2el = _relation_IgJl2Il(c2e, e2c)
    vcl2vl = _relation_IgJl2Il(c2v, v2c)
    return Mesh3Daux(f2v, c2e, c2v, v2e, v2f, v2c, e2f, e2c, f2c, vel2vl, vfl2vl, vcl2vl, efl2el, ecl2el, fcl2fl)
end
