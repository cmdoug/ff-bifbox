# moderescale.md
Author: Chris Douglas ([@cmdoug](https://github.com/cmdoug)) [christopher.douglas@duke.edu](mailto:christopher.douglas@duke.edu)

MUST BE RUN WITH 1 MPI PROCESS

```freefem
load "iovtk"
include "settings.idp"
include "macros_bifbox.idp"
// arguments
string meshin = getARGV("-mi", ""); // input meshfile
string filein = getARGV("-fi", "");
string fileout = getARGV("-fo", "");
real timeshift = getARGV("-timeshift", 0.0);
real phaseshift = getARGV("-phaseshift", 0.0);
assert(mpisize==1);
//string basefileroot, basefileext = parsefilename(basefilein, basefileroot);
string fileroot, fileext = parsefilename(filein, fileroot); //extract file name and extension
if(filein != "" && meshin == "") meshin = readmeshname(workdir + filein); // get mesh file
string meshroot, meshext = parsefilename(meshin, meshroot);
// Load mesh, make FE basis
Th = readmeshN(workdir + meshin);
Thg = Th;
DmeshCreate(Th);
restu = restrict(XMh, XMhg, n2o);
XMh defu(ub), defu(ub2);
XMh<complex> defu(um), defu(um2), defu(um3);

real omega;
if(fileext == "porb") {
    int Nh;
    complex[int, int] qh(um[].n, Nh);
    ub[] = loadporb(fileroot, meshin, qh, sym, omega, Nh);
    um[] = 0.0;
    for (int nh = 0; nh < Nh; nh++) {
        if(timeshift > 0)
            qh(:, nh) *= exp(1i*(nh+1)*omega*timeshift);
        else if(phaseshift > 0)
            qh(:, nh) *= exp(1i*(nh+1)*phaseshift);
        um[] += (2.0*qh(:, nh));
    }
}

real amp = getARGV("-amp", 1.0);
ub[] += amp*um[].re;

if(timeshift == 0 && phaseshift == 0) {
    savebase(fileout, "", meshin, true, true);
}
else {
    if(phaseshift > 0.0 && timeshift == 0) timeshift = phaseshift*omega;
    savetdns(fileout, fileout, meshin, filein, timeshift, false, true);
}
```