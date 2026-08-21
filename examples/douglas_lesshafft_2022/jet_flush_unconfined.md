# jet_flush_unconfined.md
Author: Chris Douglas ([@cmdoug](https://github.com/cmdoug)) [christopher.douglas@duke.edu](mailto:christopher.douglas@duke.edu)

This file can be used with FreeFEM to create a mesh for the flush-mounted, radially-unconfined jet as in [Douglas & Lesshafft. JFM, (2022)](https://doi.org/10.1017/jfm.2022.589).

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

string meshout = getARGV("-mo", "jet_flush_unconfined.msh"); // mesh filename
if(meshout.rfind(".msh") < 0) meshout = meshout + ".msh"; // add extension if not provided

// Define borders
//          o-4-o--..__
//          5   9      "=..
//          o-6-o          3
//              7           \
//  o-----8-----o            \
//  1                         |
//  o------------2------------o
border C01(t=0, 1){x = -l; y = 0.5*(1.0 - t); label = BCinflow;}
border C02(t=0, 1){x = -l + (l + Rinf)*t; y = 0.0; label = BCaxis;}
border C03(t=0, 1){x = Rinf*cos(pi/2.0*t); y = Rinf*sin(pi/2.0*t); label=BCopen;}
//border C04(t=0, 1){x = -Xw*t; y = Rinf; label = BCopen;}
//border C05(t=0, 1){x = -Xw; y = Rinf - (Rinf - 0.5 - deltaw)*t; label = BCwall;}
//border C06(t=0, 1){x = Xw*(t - 1.0); y = 0.5 + deltaw; label = BCwall;}
border C07(t=0, 1){x = 0.0; y = 0.5 + deltaw*(1.0 - t); label = BCwall;}
border C08(t=0, 1){x = -l*t; y = 0.5; label = BCpipe;}
border C09(t=0, 1){x = 0.0; y = Rinf - (Rinf - 0.5 - deltaw)*t; label = BCwall;}
// Assemble mesh
mesh Thg = buildmesh(C01(0.5*n0) + C02((Rinf + l)*n1) + C03(Rinf*pi/2.0*n2) 
										 + C09((Rinf - 0.5 - deltaw)*n1) + C07(deltaw*n0) + C08(l*n0));

plot(Thg,wait=1);
int[int] meshlabels = labels(Thg);
cout << "\tMesh: " << Thg.nv << " vertices, " << Thg.nt << " elements, " << Thg.nbe << " boundary elements, " << meshlabels.n << " labeled boundaries." << endl;
cout << "  Saving mesh '" + meshout + "'." << endl;
savemesh(Thg, meshout);
```