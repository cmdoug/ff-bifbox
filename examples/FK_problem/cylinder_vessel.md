# cylinder_vessel.md
Author: Mario Napieralski-Fernandez ([@marionapifm](https://github.com/marionapifm)) [mnapiera@pa.uc3m.es](mailto:mnapiera@pa.uc3m.es)

This file can be used with FreeFEM to create a mesh for the cylinder vessel for Frank-Kamenetskii problem

```freefem
assert(mpisize == 1); // Must be run with 1 processor
include "settings.idp"


int n = 350;


string meshout = getARGV("-mo","cylinder_vessel.msh"); 
if(meshout.rfind(".msh") <0) meshout = meshout + ".msh"; // add extension if not provided

// mesh filename, if not provided, defaults to "cylinder_vessel.msh"

border circle(t = 0, 2*pi) { x = cos(t); y = sin(t); label = BCwall; }
mesh Thg = buildmesh(circle(n));
int[int] meshlabels = labels(Thg);
cout << "\tMesh: " << Thg.nv << " vertices, " << Thg.nt << " elements, " << Thg.nbe << " boundary elements, " << meshlabels.n << " labeled boundaries." << endl;
cout << "  Saving mesh '" + meshout + "'." << endl;
savemesh(Thg, meshout);
plot(Thg);
```