# jet_flush_confined.md
Author: Chris Douglas ([@cmdoug](https://github.com/cmdoug)) [christopher.douglas@duke.edu](mailto:christopher.douglas@duke.edu)

This file can be used with FreeFEM to create a mesh for the flush-mounted, radially-confined jet as in [Douglas & Lesshafft. JFM, (2022)](https://doi.org/10.1017/jfm.2022.589).

```freefem
assert(mpisize == 1); // Must be run with 1 processor
include "settings.idp"
real n0 = 25;
real n1 = 5;
real n2 = 1;
real l = 4.0;
real deltaw = getARGV("-deltaw", 1.0e-5);
real Xw = 1.0;
real Rw = 8.0;
real Rinf = 40.0;
real Xinf = 100.0;

string meshout = getARGV("-mo", "jet_flush_confined.msh"); // mesh filename
if(meshout.rfind(".msh") < 0) meshout = meshout + ".msh"; // add extension if not provided

// Define borders
//          o-6-o------5------o
//          7   11            |
//          o-8-o             4
//              9             |
//  o-----10----o-----12------o
//  1                         3
//  o------------2------------o
border C01(t=0, 1){x = -l; y = 0.5*(1.0 - t); label = BCinflow;}
border C02(t=0, 1){x = -l + (l + Xinf)*t; y = 0.0; label = BCaxis;}
border C03(t=0, 1){x = Xinf; y = (0.5 + deltaw)*t; label = BCopen;}
border C04(t=0, 1){x = Xinf; y = 0.5 + deltaw + (Rw - 0.5 - deltaw)*t; label = BCopen;}
border C05(t=0, 1){x = Xinf*(1.0 - t); y = Rw; label = BCwall;}
//border C06(t=0, 1){x = -Xw*t; y = Rw; label = BCwall;}
//border C07(t=0, 1){x = -Xw; y = Rw - (Rw-0.5-deltaw)*t; label = BCwall;}
//border C08(t=0, 1){x = Xw*(t - 1.0); y = 0.5 + deltaw; label = BCwall;}
border C09(t=0, 1){x = 0.0; y = 0.5 + deltaw*(1.0 - t); label = BCwall;}
border C10(t=0, 1){x = -l*t; y = 0.5; label = BCpipe;}
border C11(t=0, 1){x = 0; y = Rw - (Rw - 0.5 - deltaw)*t; label = BCwall;}
border C12(t=0, 1){x = Xinf*t; y = 0.5+deltaw; label = BCinternal;}
// Assemble mesh
mesh Thg = buildmesh(C01(0.5*n0) + C02((Xinf + l)*n1) + C03((0.5 + deltaw)*n2) 
                     + C04((Rw - 0.5 - deltaw)*n2) + C05(Xinf*n2) + C09(deltaw*n0)
										 + C10(l*n0) + C11((Rw - 0.5 - deltaw)*n1) + C12(Xinf*n0));

plot(Thg,wait=1);
int[int] meshlabels = labels(Thg);
cout << "\tMesh: " << Thg.nv << " vertices, " << Thg.nt << " elements, " << Thg.nbe << " boundary elements, " << meshlabels.n << " labeled boundaries." << endl;
cout << "  Saving mesh '" + meshout + "'." << endl;
savemesh(Thg, meshout);
```