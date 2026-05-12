# extend_brusselator.md
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
int Nh0, Nh1;
real[int] sym0(sym.n), sym1(sym.n);
real omega;
assert(mpisize==1);

string fileroot, fileext = parsefilename(filein, fileroot); //extract file name and extension
string outfileext = parsefilename(fileout, fileout); // trim extension from output file, if given
string meshin = readmeshname(workdir + filein); // get mesh file
string meshroot, meshext = parsefilename(meshin, meshroot);
string symmetry = meshroot(meshroot.rfind("_")+1:meshroot.length-1); // get file extension
meshin = meshroot(0:meshroot.rfind("_")-1); // get file root less symmetry and file extension
meshin = meshin + "_full." + meshext; // add mesh extension for full domain

// Load mesh, make FE basis
Th = readmeshN(workdir + meshin);
Thg = Th;
DmeshCreate(Th);
restu = restrict(XMh, XMhg, n2o);
XMh defu(ub);
XMh<complex> defu(um), defu(um2), defu(um3);
complex[int, int] uh(um[].n, max(1, Nh0)), umh(um[].n, max(1, 2*Nh1));
complex eigenvalue;
if (fileext == "porb") {
  ub[] = loadporb(fileroot, meshin, uh, sym0, omega, Nh0);
}
else if (fileext == "floq"){
  um[] = loadfloq(fileroot, meshin, umh, sym1, eigenvalue, sym0, omega, Nh1);
  string porbfileroot, porbfilein = readbasename(workdir + filein);
  parsefilename(porbfilein, porbfileroot);
  ub[] = loadporb(porbfileroot, meshin, uh, sym0, omega, Nh0);
}
int isasym;
if (symmetry == "eighth" || symmetry == "Qy" || symmetry == "Qz" || symmetry == "Hyz"){
  defu(ub) = [(x >= 0)*ub(x, y, z) + (x < 0)*ub(-x, y, z), 
              (x >= 0)*ubY(x, y, z) + (x < 0)*ubY(-x, y, z)];
  defu(um) = [(x >= 0)*um(x, y, z) + (sym1(0) ? -1 : 1)*(x < 0)*um(-x, y, z), 
              (x >= 0)*umY(x, y, z) + (sym1(0) ? -1 : 1)*(x < 0)*umY(-x, y, z)];
  for (int ii = 0; ii < Nh0; ii++) {
    isasym = (sym0(0)*(ii+1) % 2);
    um2[] = uh(:, ii);
    defu(um3) = [(x >= 0)*um2(x, y, z) + (isasym ? -1 : 1)*(x < 0)*um2(-x, y, z),
                 (x >= 0)*um2Y(x, y, z) + (isasym ? -1 : 1)*(x < 0)*um2Y(-x, y, z)];
    uh(:, ii) = um3[];
  }
  for (int ii = 0; ii < Nh1; ii++) {
    isasym = (sym1(0)*(ii) % 2);
    um2[] = umh(:, 2*ii);
    defu(um3) = [(x >= 0)*um2(x, y, z) + (isasym ? -1 : 1)*(x < 0)*um2(-x, y, z),
                 (x >= 0)*um2Y(x, y, z) + (isasym ? -1 : 1)*(x < 0)*um2Y(-x, y, z)];
    umh(:, 2*ii) = um3[];
    um2[] = umh(:, 2*ii+1);
    defu(um3) = [(x >= 0)*um2(x, y, z) + (isasym ? -1 : 1)*(x < 0)*um2(-x, y, z),
                 (x >= 0)*um2Y(x, y, z) + (isasym ? -1 : 1)*(x < 0)*um2Y(-x, y, z)];
    umh(:, 2*ii+1) = um3[];
  }
}
if (symmetry == "eighth" || symmetry == "Qx" || symmetry == "Qz" || symmetry == "Hxz"){
  defu(ub) = [(y >= 0)*ub(x, y, z) + (y < 0)*ub(x, -y, z),
              (y >= 0)*ubY(x, y, z) + (y < 0)*ubY(x, -y, z)];
  defu(um) = [(y >= 0)*um(x, y, z) + (sym1(1) ? -1 : 1)*(y < 0)*um(x, -y, z), 
              (y >= 0)*umY(x, y, z) + (sym1(1) ? -1 : 1)*(y < 0)*umY(x, -y, z)];
  for (int ii = 0; ii < Nh0; ii++) {
    um2[] = uh(:, ii);
    isasym = (sym0(1)*(ii+1) % 2);
    defu(um3) = [(y >= 0)*um2(x, y, z) + (isasym ? -1 : 1)*(y < 0)*um2(x, -y, z),
                 (y >= 0)*um2Y(x, y, z) + (isasym ? -1 : 1)*(y < 0)*um2Y(x, -y, z)];
    uh(:, ii) = um3[];
  }
  for (int ii = 0; ii < Nh1; ii++) {
    isasym = (sym1(1)*(ii) % 2);
    um2[] = umh(:, 2*ii);
    defu(um3) = [(y >= 0)*um2(x, y, z) + (isasym ? -1 : 1)*(y < 0)*um2(x, -y, z),
                 (y >= 0)*um2Y(x, y, z) + (isasym ? -1 : 1)*(y < 0)*um2Y(x, -y, z)];
    umh(:, 2*ii) = um3[];
    um2[] = umh(:, 2*ii+1);
    defu(um3) = [(y >= 0)*um2(x, y, z) + (isasym ? -1 : 1)*(y < 0)*um2(x, -y, z),
                 (y >= 0)*um2Y(x, y, z) + (isasym ? -1 : 1)*(y < 0)*um2Y(x, -y, z)];
    umh(:, 2*ii+1) = um3[];
  }
}
if (symmetry == "eighth" || symmetry == "Qx" || symmetry == "Qy" || symmetry == "Hxy"){
  defu(ub) = [(z >= 0)*ub(x, y, z) + (z < 0)*ub(x, y, -z),
              (z >= 0)*ubY(x, y, z) + (z < 0)*ubY(x, y, -z)];
  defu(um) = [(z >= 0)*um(x, y, z) + (sym1(2) ? -1 : 1)*(z < 0)*um(x, y, -z), 
              (z >= 0)*umY(x, y, z) + (sym1(2) ? -1 : 1)*(z < 0)*umY(x, y, -z)];
  for (int ii = 0; ii < Nh0; ii++) {
    isasym = (sym0(2)*(ii+1) % 2);
    um2[] = uh(:, ii);
    defu(um3) = [(z >= 0)*um2(x, y, z) + (isasym ? -1 : 1)*(z < 0)*um2(x, y, -z),
                 (z >= 0)*um2Y(x, y, z) + (isasym ? -1 : 1)*(z < 0)*um2Y(x, y, -z)];
    uh(:, ii) = um3[];
  }
  for (int ii = 0; ii < Nh1; ii++) {
    isasym = (sym1(2)*(ii) % 2);
    um2[] = umh(:, 2*ii);
    defu(um3) = [(z >= 0)*um2(x, y, z) + (isasym ? -1 : 1)*(z < 0)*um2(x, y, -z),
                 (z >= 0)*um2Y(x, y, z) + (isasym ? -1 : 1)*(z < 0)*um2Y(x, y, -z)];
    umh(:, 2*ii) = um3[];
    um2[] = umh(:, 2*ii+1);
    defu(um3) = [(z >= 0)*um2(x, y, z) + (isasym ? -1 : 1)*(z < 0)*um2(x, y, -z),
                 (z >= 0)*um2Y(x, y, z) + (isasym ? -1 : 1)*(z < 0)*um2Y(x, y, -z)];
    umh(:, 2*ii+1) = um3[];
  }
}
saveporb(fileout + "start", "", meshin, sym0, omega, Nh0, true, true);
if (fileext == "floq") {
Nh0 = min(Nh0, Nh1);
ub[] += 2.0*amplitude*um[].re;
for (int ii = 0; ii < Nh0; ii++) {
 uh(:, ii) += amplitude*umh(:,2*ii);
 complex[int] qq = amplitude*umh(:,2*ii+1);
 uh(:, ii) += conj(qq);
}
saveporb(fileout + "guess", "", meshin, sym1, omega, Nh0, true, true);
ub[] = 2.0*amplitude*um[].re;
for (int ii = 0; ii < Nh0; ii++) {
 uh(:, ii) = amplitude*umh(:,2*ii);
 complex[int] qq = amplitude*umh(:,2*ii+1);
 uh(:, ii) += conj(qq);
}
saveporb(fileout + "branch", "", meshin, sym1, omega, Nh0, true, true);
}
```