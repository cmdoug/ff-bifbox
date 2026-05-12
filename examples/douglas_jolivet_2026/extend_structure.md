# extend_structure.md
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
meshin = meshroot(0:meshroot.rfind("_")-1); // get file root less symmetry and file extension
meshin = meshin + "_A." + meshext; // add mesh extension for full domain

// Load mesh, make FE basis
Th = readmeshN(workdir + meshin);
Thg = Th;
DmeshCreate(Th);
restu = restrict(XMh, XMhg, n2o);
XMh defu(ub);
XMh defu(um), defu(uma);

if (fileext == "base") ub[] = loadbase(fileroot, meshin);
else if(fileext == "fold") {
  real[string] alpha;
  real beta;
  ub[] = loadfold(fileroot, meshin, um[], uma[], alpha, beta);
}
else if (fileext == "hopf") {
  complex eigenvalue;
  real omega;
  complex[string] alpha;
  complex beta;
  complex[int] qm(um[].n), qma(um[].n);
  ub[] = loadhopf(fileroot, meshin, qm, qma, sym, omega, alpha, beta);
  um[] = 2.*qm.re;
  uma[] = 2.*qma.re;
}
if (symmetry == "S" || symmetry == "Hx"){
  defu(ub) = [(x <= 0)*ub(x, y, z) - (x > 0)*ub(-x, y, z), 
              (x <= 0)*uby(x, y, z) + (x > 0)*uby(-x, y, z),
              (x <= 0)*ubz(x, y, z) + (x > 0)*ubz(-x, y, z)];
  defu(um) = [(x <= 0)*um(x, y, z) - (sym(0) ? -1 : 1)*(x > 0)*um(-x, y, z),
              (x <= 0)*umy(x, y, z) + (sym(0) ? -1 : 1)*(x > 0)*umy(-x, y, z),
              (x <= 0)*umz(x, y, z) + (sym(0) ? -1 : 1)*(x > 0)*umz(-x, y, z)];
}

if (symmetry == "S" || symmetry == "Hy"){
  defu(ub) = [(y <= 0)*ub(x, y, z) + (y > 0)*ub(x, -y, z),
              (y <= 0)*uby(x, y, z) - (y > 0)*uby(x, -y, z),
              (y <= 0)*ubz(x, y, z) + (y > 0)*ubz(x, -y, z)];
  defu(um) = [(y <= 0)*um(x, y, z) + (sym(1) ? -1 : 1)*(y > 0)*um(x, -y, z),
              (y <= 0)*umy(x, y, z) - (sym(1) ? -1 : 1)*(y > 0)*umy(x, -y, z),
              (y <= 0)*umz(x, y, z) + (sym(1) ? -1 : 1)*(y > 0)*umz(x, -y, z)];
}

savebase(fileout + "start", "", meshin, true, true);
ub[] += amplitude*um[];
savebase(fileout + "guess", "", meshin, true, true);
ub[] = amplitude*um[];
savebase(fileout + "branch", "", meshin, true, true);
```