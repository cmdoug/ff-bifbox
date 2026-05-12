# extend_cylinder.md
Author: Chris Douglas ([@cmdoug](https://github.com/cmdoug)) [christopher.douglas@duke.edu](mailto:christopher.douglas@duke.edu)

MUST BE RUN WITH 1 MPI PROCESS

```freefem
load "iovtk"
include "settings.idp"
include "macros_bifbox.idp"
// arguments
string filein = getARGV("-fi", "");
string fileout = getARGV("-fo", "");
real amplitude = getARGV("-amp", 1.0);

assert(mpisize==1);

string fileroot, fileext = parsefilename(filein, fileroot); //extract file name and extension
string outfileext = parsefilename(fileout, fileout); // trim extension from output file, if given
if(fileext == "mode" || fileext == "resp" || fileext == "rslv" || fileext == "tdls" || fileext == "floq"){
  filein = readbasename(workdir + filein);
  fileext = parsefilename(filein, fileroot);
}
string meshin = readmeshname(workdir + filein); // get mesh file
string meshroot, meshext = parsefilename(meshin, meshroot);
string symmetry = meshroot(meshroot.rfind("_")+1:meshroot.length-1); // get file extension

// Load mesh, make FE basis
Th = readmeshN(workdir + meshin);
meshin = meshroot + "full." + meshext; // add mesh extension for full domain
mesh Th2 = movemesh(Th,[x,-y]);
Th = Th + Th2;
Th = change(Th, rmInternalEdges=1);
savemesh(Th, workdir + meshin);
Thg = Th;
DmeshCreate(Th);
restu = restrict(XMh, XMhg, n2o);
XMh defu(ub);
XMh defu(um), defu(uma);

int  Nh = getARGV("-Nh", 1);
complex[int, int] uh(ub[].n, Nh);
real omega;
ub[] = loadporb(fileroot, meshin, uh, sym, omega, Nh);

defu(ub) = [(y >= 0)*ub(x, y) + (y < 0)*ub(x, -y),
            (y >= 0)*uby(x, y) - (y < 0)*uby(x, -y),
            (y >= 0)*ubT(x, y) + (y < 0)*ubT(x, -y),
            (y >= 0)*ubp(x, y) + (y < 0)*ubp(x, -y)];
for (int ii = 0; ii < Nh; ii++) {
  uma[] = 2.*uh(:, ii).re;
  int isasym = sym(0)*(ii+1);
  defu(um) = [(y >= 0)*uma(x, y) + (isasym ? -1 : 1)*(y < 0)*uma(x, -y),
              (y >= 0)*umay(x, y) - (isasym ? -1 : 1)*(y < 0)*umay(x, -y),
              (y >= 0)*umaT(x, y) + (isasym ? -1 : 1)*(y < 0)*umaT(x, -y),
              (y >= 0)*umap(x, y) + (isasym ? -1 : 1)*(y < 0)*umap(x, -y)];
  ub[] += um[];
}
savebase(fileout, "", meshin, true, true);
```