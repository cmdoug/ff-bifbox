# vortex.md
Author: Chris Douglas ([@cmdoug](https://github.com/cmdoug)) [christopher.douglas@duke.edu](mailto:christopher.douglas@duke.edu)

This file can be used to create a mesh for the Grabowski–Berger vortex problem as in [Meliga, Gallaire, & Chomaz. JFM. (2012)](https://doi.org/10.1017/jfm.2012.93).

```freefem
assert(mpisize == 1); // Must be run with 1 processor
include "settings.idp"
real n0 = 10;
real n1 = 3;
real n2 = 1;
real xmax = 120.0;
real rmax = 60.0;

string meshout = getARGV("-mo", "vortex.msh"); // mesh filename
if(meshout.rfind(".msh") < 0) meshout = meshout + ".msh"; // add extension if not provided

// Define borders
//  o-------------3-----------o
//  |                         |
//  4                         |
//  |                         2
//  |                         |
//  o----------1--------------o
border C01(t=0, 1){x=xmax*t; y=0; label=BCaxis;}
border C02(t=0, 1){x=xmax; y=rmax*t; label=BCopen;}
border C03(t=0, 1){x=xmax*(1-t); y=rmax; label=BClat;}
border C04(t=0, 1){x=0.0; y=rmax*(1-t); label=BCinflow;}
// Assemble mesh
mesh Thg = buildmesh(C01(n1*xmax) + C02(n2*rmax) + C03(n2*xmax) + C04(n0*rmax));

int[int] meshlabels = labels(Thg);
cout << "\tMesh: " << Thg.nv << " vertices, " << Thg.nt << " elements, " << Thg.nbe << " boundary elements, " << meshlabels.n << " labeled boundaries." << endl;
cout << "  Saving mesh '" + meshout + "'." << endl;
savemesh(Thg, meshout);
```