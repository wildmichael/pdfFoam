# Auto-Ignition Burner

This is a pure demonstration show-case. The combustion model is very simple,
refer to [[Kulkarni2012]](../../doc/references.md#Kulkarni2012) for more
information. The tabulated data can be generated with the scripts
[genAutoIgnitionTable](tools/genAutoIgnitionTable) and
[convertAutoIgnitionTable](tools/convertAutoIgnitionTable). The geometry has
been created with Catia and the mesh was generated with Icem CFD. First, it was
converted to Fluent 6 format, and then the utility `fluent3DMeshToFoam` was
used to convert the mesh to OpenFOAM format. Finally, the mesh was rotated such
that the stream-wise direction aligns with the x-coordinate direction. The last
two steps are performed by the [makeMesh](makeMesh) script.
