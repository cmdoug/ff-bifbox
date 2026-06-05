# mesh_brusselator.md
Author: Chris Douglas ([@cmdoug](https://github.com/cmdoug)) [christopher.douglas@duke.edu](mailto:christopher.douglas@duke.edu)

This file is used to create initial meshes for the example of the 3-D Brusselator from [Douglas, Jolivet. CPC. (2026)](https://doi.org/10.1016/j.cpc.2026.110221).

```freefem
assert(mpisize == 1); // Must be run with 1 processor
include "settings.idp" // so that boundary labels match the convention in the settings file
int n2 = ceil(getARGV("-n", 10)/2);
string meshout = getARGV("-mo", "cube"); // mesh filename
if(meshout.rfind(".mesh") > 0) meshout = meshout(0:meshout.length - 6); // add extension if not provided

// ASSEMBLE CUBE MESH USING REFLECTIONS
// base 1/8 mesh
mesh3 Th8 = cube(n2, n2, n2, flags = 5); // create mesh
Th8 = movemesh(Th8, [x/2, y/2, z/2]);
int[int] labs = [1, BCsY, 2, BCnX, 3, BCnY, 4, BCsX, 5, BCsZ, 6, BCnZ];
Th8 = change(Th8, label = labs);
int[int] meshlabels = labels(Th8);
cout << "\tMesh: " << Th8.nv << " vertices, " << Th8.nt << " elements, " << Th8.nbe << " boundary elements, " << meshlabels.n << " labeled boundaries." << endl;
cout << "  Saving mesh '" + meshout + "_eighth.mesh'." << endl;
savemesh(Th8, meshout + "_eighth.mesh");
// 1/4 mesh -- reflected across x
mesh3 Th4x = movemesh(Th8, [-x, y, z], orientation = -1); // move mesh to center
Th4x = Th8 + Th4x;
labs = [1, BCsY, 2, BCnX, 3, BCnY, 4, BCdX, 5, BCsZ, 6, BCnZ];
Th4x = change(Th4x, label = labs, rmInternalFaces = true);
meshlabels = labels(Th4x);
cout << "\tMesh: " << Th4x.nv << " vertices, " << Th4x.nt << " elements, " << Th4x.nbe << " boundary elements, " << meshlabels.n << " labeled boundaries." << endl;
cout << "  Saving mesh '" + meshout + "_Qx.mesh'." << endl;
savemesh(Th4x, meshout + "_Qx.mesh");
// 1/4 mesh -- reflected across y 
mesh3 Th4y = movemesh(Th8, [x, -y, z], orientation = -1); // move mesh to center
Th4y = Th8 + Th4y;
labs = [1, BCdY, 2, BCnX, 3, BCnY, 4, BCsX, 5, BCsZ, 6, BCnZ];
Th4y = change(Th4y, label = labs, rmInternalFaces = true);
meshlabels = labels(Th4y);
cout << "\tMesh: " << Th4y.nv << " vertices, " << Th4y.nt << " elements, " << Th4y.nbe << " boundary elements, " << meshlabels.n << " labeled boundaries." << endl;
cout << "  Saving mesh '" + meshout + "_Qy.mesh'." << endl;
savemesh(Th4y, meshout + "_Qy.mesh");
// 1/4 mesh -- reflected across z
mesh3 Th4z = movemesh(Th8, [x, y, -z], orientation = -1); // move mesh to center
Th4z = Th8 + Th4z;
labs = [1, BCsY, 2, BCnX, 3, BCnY, 4, BCsX, 5, BCdZ, 6, BCnZ];
Th4z = change(Th4z, label = labs, rmInternalFaces = true);
meshlabels = labels(Th4z);
cout << "\tMesh: " << Th4z.nv << " vertices, " << Th4z.nt << " elements, " << Th4z.nbe << " boundary elements, " << meshlabels.n << " labeled boundaries." << endl;
cout << "  Saving mesh '" + meshout + "_Qz.mesh'." << endl;
savemesh(Th4z, meshout + "_Qz.mesh");
// 1/2 mesh -- reflected across x,y
mesh3 Th2xy = movemesh(Th4x, [x, -y, z], orientation = -1); // move mesh to center
Th2xy = Th4x + Th2xy;
labs = [1, BCdY, 2, BCnX, 3, BCnY, 4, BCdX, 5, BCsZ, 6, BCnZ];
Th2xy = change(Th2xy, label = labs, rmInternalFaces = true);
meshlabels = labels(Th2xy);
cout << "\tMesh: " << Th2xy.nv << " vertices, " << Th2xy.nt << " elements, " << Th2xy.nbe << " boundary elements, " << meshlabels.n << " labeled boundaries." << endl;
cout << "  Saving mesh '" + meshout + "_Hxy.mesh'." << endl;
savemesh(Th2xy, meshout + "_Hxy.mesh");
// 1/2 mesh -- reflected across x,z
mesh3 Th2xz = movemesh(Th4x, [x, y, -z], orientation = -1); // move mesh to center
Th2xz = Th4x + Th2xz;
labs = [1, BCsY, 2, BCnX, 3, BCnY, 4, BCdX, 5, BCdZ, 6, BCnZ];
Th2xz = change(Th2xz, label = labs, rmInternalFaces = true);
meshlabels = labels(Th2xz);
cout << "\tMesh: " << Th2xz.nv << " vertices, " << Th2xz.nt << " elements, " << Th2xz.nbe << " boundary elements, " << meshlabels.n << " labeled boundaries." << endl;
cout << "  Saving mesh '" + meshout + "_Hxz.mesh'." << endl;
savemesh(Th2xz, meshout + "_Hxz.mesh");
// 1/2 mesh -- reflected across y,z
mesh3 Th2yz = movemesh(Th4y, [x, y, -z], orientation = -1); // move mesh to center
Th2yz = Th4y + Th2yz;
labs = [1, BCdY, 2, BCnX, 3, BCnY, 4, BCsX, 5, BCdZ, 6, BCnZ];
Th2yz = change(Th2yz, label = labs, rmInternalFaces = true);
meshlabels = labels(Th2yz);
cout << "\tMesh: " << Th2yz.nv << " vertices, " << Th2yz.nt << " elements, " << Th2yz.nbe << " boundary elements, " << meshlabels.n << " labeled boundaries." << endl;
cout << "  Saving mesh '" + meshout + "_Hyz.mesh'." << endl;
savemesh(Th2yz, meshout + "_Hyz.mesh");
// full mesh -- reflected across x,y,z
mesh3 Thg = movemesh(Th2xy, [x, y, -z], orientation = -1); // move mesh to center
Thg = Thg + Th2xy;
labs = [1, BCdY, 2, BCnX, 3, BCnY, 4, BCdX, 5, BCdZ, 6, BCnZ];
Thg = change(Thg, label = labs, rmInternalFaces = true);
meshlabels = labels(Thg);
cout << "\tMesh: " << Thg.nv << " vertices, " << Thg.nt << " elements, " << Thg.nbe << " boundary elements, " << meshlabels.n << " labeled boundaries." << endl;
cout << "  Saving mesh '" + meshout + "_full.mesh'." << endl;
savemesh(Thg, meshout + "_full.mesh");
```