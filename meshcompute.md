# meshcompute.md
Author: Chris Douglas ([@cmdoug](https://github.com/cmdoug)) [christopher.douglas@duke.edu](mailto:christopher.douglas@duke.edu)

## EXAMPLE USAGE:
### Initialize with solution from one file
```sh
ff-mpirun -np 1 meshcompute.md -fi <FILEIN> -mo <MESHOUT>
```

### Initialize with solutions from two files
```sh
ff-mpirun -np 1 meshcompute.md -fi <FILEIN> -fi2 <FILEIN2> -mo <FILEOUT>
```

NOTE: This file should not be changed unless you know what you're doing.

SEE ALSO: [basecompute.md](./basecompute.md), [modecompute.md](./modecompute.md), [foldcompute.md](./foldcompute.md), [hopfcompute.md](./hopfcompute.md), [porbcompute.md](./porbcompute.md)

```freefem
include "settings.idp"
include "macros_bifbox.idp"
// arguments
string meshin = getARGV("-mi", ""); // input meshfile with extension
string meshout = getARGV("-mo", ""); // output meshfile without extension
string filein = getARGV("-fi", ""); // input file with extension
string filein2 = getARGV("-fi2", ""); // input file with extension
string filein3 = getARGV("-fi3", ""); // input file with extension
string filein4 = getARGV("-fi4", ""); // input file with extension
string filein5 = getARGV("-fi5", ""); // input file with extension
string filein6 = getARGV("-fi6", ""); // input file with extension
string filein7 = getARGV("-fi7", ""); // input file with extension
string filein8 = getARGV("-fi8", ""); // input file with extension
string filein9 = getARGV("-fi9", ""); // input file with extension
string filein0 = getARGV("-fi0", ""); // input file with extension
string filein1 = getARGV("-fi1", ""); // input file with extension
string adaptto = getARGV("-adaptto", "b");

int filecount = (filein!="") + (filein2!="") + (filein3!="") + (filein4!="") + (filein5!="") + (filein6!="") + (filein7!="") + (filein8!="") + (filein9!="") + (filein0!="") + (filein1!="");
assert(filecount>0);
string[int] fileins(filecount), fileroots(filecount), fileexts(filecount);
int ii = 0;
if(filein!="") { fileins[ii] = filein; ii++; }
if(filein2!="") { fileins[ii] = filein2; ii++; }
if(filein3!="") { fileins[ii] = filein3; ii++; }
if(filein4!="") { fileins[ii] = filein4; ii++; }
if(filein5!="") { fileins[ii] = filein5; ii++; }
if(filein6!="") { fileins[ii] = filein6; ii++; }
if(filein7!="") { fileins[ii] = filein7; ii++; }
if(filein8!="") { fileins[ii] = filein8; ii++; }
if(filein9!="") { fileins[ii] = filein9; ii++; }
if(filein0!="") { fileins[ii] = filein0; ii++; }
if(filein1!="") { fileins[ii] = filein1; ii++; }

int nuvec = 0, nfvec = 0;
for (ii = 0; ii < filecount; ii++) {
  fileexts[ii] = parsefilename(fileins[ii], fileroots[ii]); //extract file name and extension
  if(meshin == "") meshin = readmeshname(workdir + fileins[ii]); // get mesh file
  if( fileexts[ii] == "base" || fileexts[ii] == "tdns" ) nuvec++;
  else if ( fileexts[ii] == "mode" || fileexts[ii] == "tdls" || fileexts[ii] == "resp" ) nuvec += 2;
  else if ( fileexts[ii] == "fold" ) nuvec += 3;
  else if ( fileexts[ii] == "hopf" || fileexts[ii] == "porb" ) nuvec += 5;
  else if ( fileexts[ii] == "foho" ) nuvec += 7;
  else if ( fileexts[ii] == "hoho" ) nuvec += 9;
  else if ( fileexts[ii] == "rslv" ) {nuvec += 2; nfvec += 2;}
}
string meshroot, meshext = parsefilename(meshin, meshroot);
parsefilename(meshout, meshout); // trim extension from output mesh, if given
meshout = meshout + "." + meshext;

// Load mesh
if (mpirank==0){ // Perform mesh adaptation (serially) on processor 0
  Th = readmeshN(workdir + meshin);
  Thg = Th;
  restu = 0:XMhg.ndof-1;
  restf = 0:Xhg.ndof-1;
  // Make finite element basis
  real[int, int] uvecs(nuvec, XMhg.ndof);
  real[int, int] fvecs(nfvec, Xhg.ndof);
  int jj = 0, kk = 0;
  for (ii = 0; ii < filecount; ii++) {
    if(fileexts[ii] == "base") {
      uvecs(jj++, :) = loadbase(fileroots[ii], meshin);
    }
    else if(fileexts[ii] == "tdns") {
      real time;
      uvecs(jj++, :) = loadtdns(fileroots[ii], meshin, time);
    }
    else if(fileexts[ii] == "mode") {
      complex eigenvalue;
      complex[int] qmg(XMhg.ndof);
      qmg = loadmode(fileroots[ii], meshin, sym, eigenvalue);
      uvecs(jj++, :) = qmg.re;
      uvecs(jj++, :) = qmg.im;
    }
    else if(fileexts[ii] == "resp") {
      real omega;
      complex[int] qmg(XMhg.ndof);
      qmg = loadresp(fileroots[ii], meshin, sym, omega);
      uvecs(jj++, :) = qmg.re;
      uvecs(jj++, :) = qmg.im;
    }
    else if(fileexts[ii] == "tdls") {
      real time;
      complex[int] qmg(XMhg.ndof);
      qmg = loadtdls(fileroots[ii], meshin, sym, time);
      uvecs(jj++, :) = qmg.re;
      uvecs(jj++, :) = qmg.im;
    }
    else if(fileexts[ii] == "fold") {
      real[string] alpha;
      real beta;
      real[int] qm(XMhg.ndof), qma(XMhg.ndof);
      uvecs(jj++, :) = loadfold(fileroots[ii], meshin, qm, qma, alpha, beta);
      if(adaptto == "bd" || adaptto == "bda") uvecs(jj++, :) = qm;
      if(adaptto == "ba" || adaptto == "bda") uvecs(jj++, :) = qma;
    }
    else if(fileexts[ii] == "hopf") {
      real omega;
      complex[string] alpha;
      complex beta;
      complex[int] qm(XMhg.ndof), qma(XMhg.ndof);
      uvecs(jj++, :) = loadhopf(fileroots[ii], meshin, qm, qma, sym, omega, alpha, beta);
      if(adaptto == "bd" || adaptto == "bda") { 
        uvecs(jj++, :) = qm.re; 
        uvecs(jj++, :) = qm.im;
      }
      if(adaptto == "ba" || adaptto == "bda") {
        uvecs(jj++, :) = qma.re;
        uvecs(jj++, :) = qma.im;
      }
    }
    else if(fileexts[ii] == "porb") {
      int Nh=1;
      real omega;
      complex[int, int] qh(XMhg.ndof, Nh);
      uvecs(jj++, :) = loadporb(fileroots[ii], meshin, qh, sym, omega, Nh);
      uvecs(jj++, :) = qh(:, 0).re;
      uvecs(jj++, :) = qh(:, 0).im;
      if (Nh > 1){
        uvecs(jj++, :) = qh(:, 1).re;
        uvecs(jj++, :) = qh(:, 1).im;
      }
    }
    else if(fileexts[ii] == "foho") {
      real omega;
      complex[string] alpha1;
      real[string] alpha2;
      complex beta1, gamma12, gamma13;
      real beta22, beta23, gamma22, gamma23;
      complex[int] q1m(XMhg.ndof), q1ma(XMhg.ndof);
      real[int] q2m(XMhg.ndof), q2ma(XMhg.ndof);
      uvecs(jj++, :) = loadfoho(fileroots[ii], meshin, q1m, q1ma, q2m, q2ma, sym, omega, alpha1, alpha2, beta1, beta22, beta23, gamma12, gamma13, gamma22, gamma23);
      if(adaptto == "bd" || adaptto == "bda") {
        uvecs(jj++, :) = q1m.re;
        uvecs(jj++, :) = q1m.im;
        uvecs(jj++, :) = q2m;
      }
      if(adaptto == "ba" || adaptto == "bda") {
        uvecs(jj++, :) = q1ma.re;
        uvecs(jj++, :) = q1ma.im;
        uvecs(jj++, :) = q2ma;
      }
    }
    else if(fileexts[ii] == "hoho") {
      real[int] sym1(sym.n), sym2(sym.n);
      real omega1, omega2;
      complex[string] alpha1, alpha2;
      complex beta1, beta2, gamma11, gamma12, gamma13, gamma21, gamma22, gamma23;
      complex[int] q1m(XMhg.ndof), q1ma(XMhg.ndof), q2m(XMhg.ndof), q2ma(XMhg.ndof);
      uvecs(jj++, :) = loadhoho(fileroots[ii], meshin, q1m, q1ma, q2m, q2ma, sym1, sym2, omega1, omega2, alpha1, alpha2, beta1, beta2, gamma11, gamma12, gamma13, gamma21, gamma22, gamma23);
      if(adaptto == "bd" || adaptto == "bda") {
        uvecs(jj++, :) = q1m.re;
        uvecs(jj++, :) = q1m.im;
        uvecs(jj++, :) = q2m.re;
        uvecs(jj++, :) = q2m.im;
      }
      if(adaptto == "ba" || adaptto == "bda") {
        uvecs(jj++, :) = q1ma.re;
        uvecs(jj++, :) = q1ma.im;
        uvecs(jj++, :) = q2ma.re;
        uvecs(jj++, :) = q2ma.im;
      }
    }
    else if(fileexts[ii] == "rslv"){
      real omega;
      real gain;
      complex[int] qmg, fmg;
      qmg = loadrslv(fileroots[ii], meshin, fmg, sym, omega, gain);
      uvecs(jj++, :) = qmg.re;
      uvecs(jj++, :) = qmg.im;
//      fvecs(kk++, :) = fmg.re;
//      fvecs(kk++, :) = fmg.im;
    }
  }
  XMhg defu(uG), defu(uG2), defu(uG3), defu(uG4), defu(uG5), defu(uG6), defu(uG7), defu(uG8), defu(uG9), defu(uG0), defu(uG1); // create private global FE functions
  uG[] = uvecs(0, :);
  if(jj > 1) uG2[] = uvecs(1, :);
  if(jj > 2) uG3[] = uvecs(2, :);
  if(jj > 3) uG4[] = uvecs(3, :);
  if(jj > 4) uG5[] = uvecs(4, :);
  if(jj > 5) uG6[] = uvecs(5, :);
  if(jj > 6) uG7[] = uvecs(6, :);
  if(jj > 7) uG8[] = uvecs(7, :);
  if(jj > 8) uG9[] = uvecs(8, :);
  if(jj > 9) uG0[] = uvecs(9, :);
  if(jj > 10) uG0[] = uvecs(10, :);
  // MESH ADAPTATION
  IFMACRO(dimension,2)
  if (jj == 1) Thg = adaptmesh(Thg, adaptu(uG), adaptmeshoptions);
  else if (jj == 2) Thg = adaptmesh(Thg, adaptu(uG), adaptu(uG2), adaptmeshoptions);
  else if (jj == 3) Thg = adaptmesh(Thg, adaptu(uG), adaptu(uG2), adaptu(uG3), adaptmeshoptions);
  else if (jj == 4) Thg = adaptmesh(Thg, adaptu(uG), adaptu(uG2), adaptu(uG3), adaptu(uG4), adaptmeshoptions);
  else if (jj == 5) Thg = adaptmesh(Thg, adaptu(uG), adaptu(uG2), adaptu(uG3), adaptu(uG4), adaptu(uG5), adaptmeshoptions);
  else if (jj == 6) Thg = adaptmesh(Thg, adaptu(uG), adaptu(uG2), adaptu(uG3), adaptu(uG4), adaptu(uG5), adaptu(uG6), adaptmeshoptions);
  else if (jj == 7) Thg = adaptmesh(Thg, adaptu(uG), adaptu(uG2), adaptu(uG3), adaptu(uG4), adaptu(uG5), adaptu(uG6), adaptu(uG7), adaptmeshoptions);
  else if (jj == 8) Thg = adaptmesh(Thg, adaptu(uG), adaptu(uG2), adaptu(uG3), adaptu(uG4), adaptu(uG5), adaptu(uG6), adaptu(uG7), adaptu(uG8), adaptmeshoptions);
  else if (jj == 9) Thg = adaptmesh(Thg, adaptu(uG), adaptu(uG2), adaptu(uG3), adaptu(uG4), adaptu(uG5), adaptu(uG6), adaptu(uG7), adaptu(uG8), adaptu(uG9), adaptmeshoptions);
  else if (jj == 10) Thg = adaptmesh(Thg, adaptu(uG), adaptu(uG2), adaptu(uG3), adaptu(uG4), adaptu(uG5), adaptu(uG6), adaptu(uG7), adaptu(uG8), adaptu(uG9), adaptu(uG0), adaptmeshoptions);
  else Thg = adaptmesh(Thg, adaptu(uG), adaptu(uG2), adaptu(uG3), adaptu(uG4), adaptu(uG5), adaptu(uG6), adaptu(uG7), adaptu(uG8), adaptu(uG9), adaptu(uG0), adaptu(uG1), adaptmeshoptions);
  ENDIFMACRO
  IFMACRO(dimension,3)
  //NOTE: 3D mesh adaptation is still under development.
  load "mshmet"
  load "mmg"
  real anisomax = getARGV("-anisomax",1.0);
  real[int] met;
  if (jj == 1) met = mshmet(Thg, adaptu(uG), normalization = getARGV("-normalization",1), aniso = bool(anisomax > 1.0), hmin = getARGV("-hmin", 1.0e-6), hmax = getARGV("-hmax", 1.0e+2), err = getARGV("-err", 1.0e-2));
  else if (jj == 2) met = mshmet(Thg, adaptu(uG), adaptu(uG2), normalization = getARGV("-normalization",1), aniso = bool(anisomax > 1.0), hmin = getARGV("-hmin", 1.0e-6), hmax = getARGV("-hmax", 1.0e+2), err = getARGV("-err", 1.0e-2));
  else if (jj == 3) met = mshmet(Thg, adaptu(uG), adaptu(uG2), adaptu(uG3), normalization = getARGV("-normalization",1), aniso = bool(anisomax > 1.0), hmin = getARGV("-hmin", 1.0e-6), hmax = getARGV("-hmax", 1.0e+2), err = getARGV("-err", 1.0e-2));
  else if (jj == 4) met = mshmet(Thg, adaptu(uG), adaptu(uG2), adaptu(uG3), adaptu(uG4), normalization = getARGV("-normalization",1), aniso = bool(anisomax > 1.0), hmin = getARGV("-hmin", 1.0e-6), hmax = getARGV("-hmax", 1.0e+2), err = getARGV("-err", 1.0e-2));
  else if (jj == 5) met = mshmet(Thg, adaptu(uG), adaptu(uG2), adaptu(uG3), adaptu(uG4), adaptu(uG5), normalization = getARGV("-normalization",1), aniso = bool(anisomax > 1.0), hmin = getARGV("-hmin", 1.0e-6), hmax = getARGV("-hmax", 1.0e+2), err = getARGV("-err", 1.0e-2));
  else if (jj == 6) met = mshmet(Thg, adaptu(uG), adaptu(uG2), adaptu(uG3), adaptu(uG4), adaptu(uG5), adaptu(uG6), normalization = getARGV("-normalization",1), aniso = bool(anisomax > 1.0), hmin = getARGV("-hmin", 1.0e-6), hmax = getARGV("-hmax", 1.0e+2), err = getARGV("-err", 1.0e-2));
  else if (jj == 7) met = mshmet(Thg, adaptu(uG), adaptu(uG2), adaptu(uG3), adaptu(uG4), adaptu(uG5), adaptu(uG6), adaptu(uG7), normalization = getARGV("-normalization",1), aniso = bool(anisomax > 1.0), hmin = getARGV("-hmin", 1.0e-6), hmax = getARGV("-hmax", 1.0e+2), err = getARGV("-err", 1.0e-2));
  else if (jj == 8) met = mshmet(Thg, adaptu(uG), adaptu(uG2), adaptu(uG3), adaptu(uG4), adaptu(uG5), adaptu(uG6), adaptu(uG7), adaptu(uG8), normalization = getARGV("-normalization",1), aniso = bool(anisomax > 1.0), hmin = getARGV("-hmin", 1.0e-6), hmax = getARGV("-hmax", 1.0e+2), err = getARGV("-err", 1.0e-2));
  else if (jj == 9) met = mshmet(Thg, adaptu(uG), adaptu(uG2), adaptu(uG3), adaptu(uG4), adaptu(uG5), adaptu(uG6), adaptu(uG7), adaptu(uG8), adaptu(uG9), normalization = getARGV("-normalization",1), aniso = bool(anisomax > 1.0), hmin = getARGV("-hmin", 1.0e-6), hmax = getARGV("-hmax", 1.0e+2), err = getARGV("-err", 1.0e-2));
  else if (jj == 10) met = mshmet(Thg, adaptu(uG), adaptu(uG2), adaptu(uG3), adaptu(uG4), adaptu(uG5), adaptu(uG6), adaptu(uG7), adaptu(uG8), adaptu(uG9), adaptu(uG0), normalization = getARGV("-normalization",1), aniso = bool(anisomax > 1.0), hmin = getARGV("-hmin", 1.0e-6), hmax = getARGV("-hmax", 1.0e+2), err = getARGV("-err", 1.0e-2));
  else met = mshmet(Thg, adaptu(uG), adaptu(uG2), adaptu(uG3), adaptu(uG4), adaptu(uG5), adaptu(uG6), adaptu(uG7), adaptu(uG8), adaptu(uG9), adaptu(uG0), adaptu(uG1), normalization = getARGV("-normalization",1), aniso = bool(anisomax > 1.0), hmin = getARGV("-hmin", 1.0e-6), hmax = getARGV("-hmax", 1.0e+2), err = getARGV("-err", 1.0e-2));
  if(anisomax > 1.0) {
    load "aniso"
    boundaniso(6, met, anisomax);
  }
  Thg = mmg3d(Thg, metric = met, hmin = getARGV("-hmin", 1.0e-6), hmax = getARGV("-hmax", 1.0e+2), hgrad = -1, verbose = verbosity-(verbosity==0));
  ENDIFMACRO
  cout << "  Saving adapted mesh '" + meshout + "' in '" + workdir + "'." << endl;
  savemesh(Thg, workdir + meshout);
}
```