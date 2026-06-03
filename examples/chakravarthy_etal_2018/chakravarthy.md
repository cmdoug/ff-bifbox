```freefem
//
// chakravarthy.geo
// Gretchen Schulke
// gretchen.schulke@duke.edu
//
// This file can be used with FreeFem to create a mesh
//

assert(mpisize == 1); // Must be run with 1 processor
include "settings.idp"
real n0 = getARGV("-n0", 5);
real n1 = getARGV("-n1", 2);
real R = getARGV("-R", 30.0);
real L = getARGV("-L", 80.0);

string meshout = getARGV("-mo", "chakravarthy.msh"); // mesh filename
if(meshout.rfind(".msh") < 0) meshout = meshout + ".msh"; // add extension if not provided

//  Define borders 
//  o-----------------3-----------------o
//  |            						|
//  4            						2
//  |            						|
//  o-----------------1-----------------o


border C01(t=0, 1){ x=L*t^2; y=0;     label=BCaxis;   }
border C02(t=0, 1){ x=L;     y=R*t;   label=BCoutlet; }
border C03(t=0, 1){ x=L*t;   y=R;     label=BCwall;   }
border C04(t=0, 1){ x=0;     y=R*t^2; label=BCinflow; }

// Assemble mesh
mesh Thg = buildmesh(C01(L*n0) + C02(R*n1) + C03(-L*n1) + C04(-R*n0));

int[int] meshlabels = labels(Thg);
cout << "\tMesh: " << Thg.nv << " vertices, " << Thg.nt << " triangles, " << Thg.nbe << " boundary elements, " << meshlabels.n << " labeled boundaries." << endl;
cout << "  Saving mesh '" + meshout + "'." << endl;
plot(Thg);
savemesh(Thg, meshout);
```