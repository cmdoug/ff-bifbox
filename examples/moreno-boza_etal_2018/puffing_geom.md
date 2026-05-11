# puffing_geom.md
Author: Daniel Moreno ([@dmorenobz](https://github.com/dmorenobz)) [damoreno@ing.uc3m.es](mailto:damoreno@ing.uc3m.es)

```freefem
assert(mpisize == 1); // Must be run with 1 processor
include "settings.idp"
real n0 = 20;
real n1 = 8;
real n2 = 1;
real xmax = 50.0;
real rmax = 12.0;
real rpipe = 1.0;

string meshout = getARGV("-mo", "jet.msh"); // mesh filename
if(meshout.rfind(".msh") < 0) meshout = meshout + ".msh"; // add extension if not provided

// Define borders
//  o------------------4-------------------o
//  |                                      |
//  5                                      |
//  |                                      3
//  o                                      |
//  1                                      |
//  o------------------2-------------------o
border C01(t=0, 1){x = 0;          y = rpipe*(1-t); 		label = BCinflow;}
border C02(t=0, 1){x = xmax*t;     y = 0; 					label = BCaxis;}
border C03(t=0, 1){x = xmax;       y = rmax*t; 				label = BCout;}
border C04(t=0, 1){x = xmax*(1-t); y = rmax; 				label = BClat;}
border C05(t=0, 1){x = 0;          y = rmax-(rmax-rpipe)*t; label = BCwall;}

// Assemble mesh
mesh Thg = buildmesh(C01(n0*rpipe) + C02(n1*xmax) + C03(n2*rmax) + C04(n2*xmax) + C05(n1*(rmax-rpipe)));

int[int] meshlabels = labels(Thg);
cout << "\tMesh: " << Thg.nv << " vertices, " << Thg.nt << " elements, " << Thg.nbe << " boundary elements, " << meshlabels.n << " labeled boundaries." << endl;
cout << "  Saving mesh '" + meshout + "'." << endl;
savemesh(Thg, meshout);
```