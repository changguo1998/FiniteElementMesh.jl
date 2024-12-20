## Naming Style

### name of mesh element

A mesh element is named as follows:

```
[(mesh tag)](element type)(element scope)
```

where `(mesh tag)` is `p`(partition mesh) or `t`(total mesh), this prefix only used when two meshs show up in same function;
`(element type)` is one of `vertex`, `edge`, `face` and `cell`;
and `(element scope)` is `l`(for local index) or `g`(for global index).
By default, the `g` can be omited.

For example, global index of edge in partition mesh can be named as `pedgeg`, and simplified as
`edge`, `eg` or `e`.

### relation matrix

matrix defined relation between two kindes of graph elements is named as follows:

```text
[(mesh tag)_](column mesh element)(row mesh element)2(matrix mesh element)
```

For example, `edgegvertl2vertg` means a matrix whose columns represent edges,
the column number is global index of edge. Each column stores global index of vertex
in this edge. The row index means the index of vertex for a single edge.

To simplify, format `xgyl2yg` can be written as `x2y`.

### Abbreviations

| Term       | Abbr. 1 char | Abbr 4 char | Expl.                                                |
| :--------- | :----------: | :---------: | :--------------------------------------------------- |
| coordinate |      x       |    coor     | 0D, a point                                          |
| vertex     |      v       |    vert     | 0D, a point                                          |
| edge       |      e       |    edge     | 1D object                                            |
| face       |      f       |    face     | 2D object                                            |
| cell       |      c       |    cell     | 3D object                                            |
| local      |      l       |    locl     | In the relation map of mesh, the index of map in row |
| global     |      g       |    glob     | Index number saved in relation map                   |
| part       |      p       |    part     | To distinguish a mesh included in another msh        |
| total      |      t       |    totl     | To distinguish a mesh containing another msh         |
