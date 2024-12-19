using FiniteElementMesh

# 1D
m1d = Mesh1D(range(0.0, 3.0, 4), [1 2 3; 2 3 4])
data1d = Mesh1Ddata((:stat, :vertex), reshape(collect(1:2*4), 2, 4))
a1d = auxillary_relation(m1d)
m1d_v = submesh(m1d, :vertex, [1, 2, 3])
m1d_e = submesh(m1d, edge=[2, 3])
data1d_e = FiniteElementMesh.submeshdata(m1d_e, m1d, data1d)

# 2D
m2d = Mesh2D(
    [
        1.0 2.0 3.0 1.0 2.0 3.0;
        0.0 0.0 0.0 1.0 1.0 1.0
    ],
    [
        1 2 1 2 3 4 5;
        2 3 4 5 6 5 6
    ],
    [
        1 2;
        4 5;
        6 7;
        3 4;
    ]
)

m2d_v = submesh(m2d, :vertex, [1, 2, 4, 5])
m2d_e = submesh(m2d, edge=[1, 3, 4, 6])
m2d_f = submesh(m2d, face=[2])

# 3D

vertex = Float64[
    0  1  0  0  0;
    0  0  1  0  0;
    0  0  0  1 -1.;
]

edge2vert = [
    1 2;
    1 3;
    2 3;
    1 4;
    2 4;
    3 4;
    1 5;
    2 5;
    3 5;
] |> permutedims

face2edge = [
    1 2 3;
    1 4 5;
    2 4 6;
    3 5 6;
    1 7 8;
    2 7 9;
    3 8 9;
] |> permutedims

cell2face = [
    1 2 3 4;
    1 5 6 7;
] |> permutedims
m3d = Mesh3D(vertex, edge2vert, face2edge, cell2face)
a3d = auxillary_relation(m3d)

FiniteElementMesh._polycoef_2d_combine_triangle(FiniteElementMesh._polycoef_1d_linear, 3)
