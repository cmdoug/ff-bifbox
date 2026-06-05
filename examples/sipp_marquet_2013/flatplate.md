# flatplate.md
Author: Daniel Moreno ([@dmorenobz](https://github.com/dmorenobz)) [damoreno@ing.uc3m.es](mailto:damoreno@ing.uc3m.es)

This file can be used with FreeFEM to create a mesh for the flat plate boundary layer as in [Sipp & Marquet. TCFD. (2013)](https://doi.org/10.1007/s00162-012-0265-y).

```freefem
assert(mpisize == 1); // Must be run with 1 processor
include "settings.idp"
real n0 = 120;
real n1 = 80;
real n2 = 50;
real xmin = -0.5;
real xmax = 1.25;
real ymax = 1.0;

string meshout = getARGV("-mo", "flatplate.msh"); // mesh filename
if(meshout.rfind(".msh") < 0) meshout = meshout + ".msh"; // add extension if not provided

//  o--------5-------------o
//  |                      |
//  1                      4
//  |                      |
//  o---2----o-----3-------o

border C01(t=0, 1){x=xmin; y=ymax*(1-t); label=BCinflow;}
border C02(t=0, 1){x=xmin*(1-t); y=0.0; label=BCslip;}
border C03(t=0, 1){x=xmax*t; y=0.0; label=BCwall;}
border C04(t=0, 1){x=xmax; y=ymax*t; label=BCopen;}
border C05(t=0, 1){x=xmax-(xmax-xmin)*t; y=ymax; label=BCslip;}

// Assemble mesh
mesh Thg = buildmesh(C01(n2*ymax) + C02(-n2*xmin) + C03(n0*xmax) + C04(n2*ymax) + C05(n1*(xmax-xmin)));

int[int] meshlabels = labels(Thg);
cout << "\tMesh: " << Thg.nv << " vertices, " << Thg.nt << " elements, " << Thg.nbe << " boundary elements, " << meshlabels.n << " labeled boundaries." << endl;
cout << "  Saving mesh '" + meshout + "'." << endl;
savemesh(Thg, meshout);
```