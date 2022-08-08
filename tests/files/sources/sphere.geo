
// Create a very coarse sphere mesh
SetFactory("OpenCASCADE");
Sphere (1) = {0.0, 0.0, 0.0, 1.0};

// Options for Saving meshes
Mesh.Format = 1; // msh
Mesh.MshFileVersion = 2.2;
Mesh.Binary = 0;

// Very coarse mesh
MeshSize{1,2} = 2.0;

Mesh 3;
SetOrder 2;

Save "sphere-quadratic-N1.msh";

RefineMesh;
SetOrder 2;
Save "sphere-quadratic-N2.msh";

RefineMesh;
SetOrder 2;
Save "sphere-quadratic-N3.msh";

RefineMesh;
SetOrder 2;
Save "sphere-quadratic-N4.msh";

