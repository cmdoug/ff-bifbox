# nozzle.md
Author: Chris Douglas ([@cmdoug](https://github.com/cmdoug)) [christopher.douglas@duke.edu](mailto:christopher.douglas@duke.edu)

This file can be used with FreeFEM to create a mesh for the turbulent swirling jet as in [Chevalier et al, JFM, (2024)](https://doi.org/10.1007/s00162-024-00704-2).

```freefem
assert(mpisize == 1); // Must be run with 1 processor
include "settings.idp"
real n0 = 25;
real n1 = 5;
real n2 = 1;
real L = getARGV("-L", 60.0);
real H = getARGV("-H", 20.0);
real w = 0.05;
real h = 1.0e-4;
real hh = 1.0e-5;

string meshout = getARGV("-mo", "nozzle.msh"); // mesh filename
if(meshout.rfind(".msh") < 0) meshout = meshout + ".msh"; // add extension if not provided

// Define borders
//  o------------4------------o
//  |                         |
//  5                         |
//  |                         |
//  o----6----o               3
//             7              |
//  o-----8-----o             |
//  1                         |
//  o------------2------------o
border C01(t=0, 1){x = -1.0; y = 1.0 - t; label = BCin1;}
border C02(t=0, 1){x = -1.0 + L*t; y = 0.0; label = BCaxis;}
border C03(t=0, 1){x = L - 1.0; y = H*t; label = BCopen;}
border C04(t=0, 1){x = -1.0 + L*(1.0 - t); y = H; label = BClat;}
border C05(t=0, 1){x = -1.0; y = H - (H - h - 1.0)*t; label = BCin2;}
border C06(t=0, 1){x = -1.0 + (1.0 - w)*t; y = 1.0 + h; label = BCwall;}
border C07(t=0, 1){x = -w*(1.0 - t); y = 1.0 + h - (h - hh)*t; label = BCwall;}
border C07a(t=0, 1){x = 0.0; y = 1.0 + hh*(1.0 - t); label = BCwall;}
border C08(t=0, 1){x = -t; y = 1.0; label = BCwall;}
// Assemble mesh
mesh Thg = buildmesh(C01(n0) + C02(L*n1) + C03(H*n2) 
                     + C04(L*n2) + C05((H - h - 1.0)*n1) + C06((1.0 - w)*n0) 
										 + C07(sqrt(w^2 + (h - hh)^2)*n0) + C07a(hh*n0) + C08(n0));

plot(Thg,wait=1);
int[int] meshlabels = labels(Thg);
cout << "\tMesh: " << Thg.nv << " vertices, " << Thg.nt << " elements, " << Thg.nbe << " boundary elements, " << meshlabels.n << " labeled boundaries." << endl;
cout << "  Saving mesh '" + meshout + "'." << endl;
savemesh(Thg, meshout);
```