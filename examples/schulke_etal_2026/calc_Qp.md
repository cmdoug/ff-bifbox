# calc_Qp.md
Author: Chris Douglas ([@cmdoug](https://github.com/cmdoug)) [christopher.douglas@duke.edu](mailto:christopher.douglas@duke.edu)

```freefem
load "iovtk"
include "settings.idp"
include "macros_bifbox.idp"

string meshin = getARGV("-mi", ""); // input meshfile with extension
string meshout = getARGV("-mo", ""); // output mesh without extension
string filein = getARGV("-fi", ""); // input file with extension
string fileout = getARGV("-fo", ""); // output file without extension
int Nh = 0;
real[int] sym0(sym.n);
real omega;


string fileroot, fileext = parsefilename(filein, fileroot); //extract file name and extension
parsefilename(fileout, fileout); // trim extension from output file, if given
if(meshin == "") meshin = readmeshname(workdir + filein); // get mesh file
string meshroot, meshext = parsefilename(meshin, meshroot);
parsefilename(meshout, meshout); // trim extension from output mesh, if given
Th = readmeshN(workdir + meshin);
Thg = Th;
DmeshCreate(Th);
restu = restrict(XMh, XMhg, n2o);
XMh defu(ub);
XMh<complex> defu(um), defu(um2), defu(um3);
complex[int, int] uh(um[].n, max(1, Nh));
assert(fileext == "porb");
ub[] = loadporb(fileroot, meshin, uh, sym0, omega, Nh);

complex[int] ik(sym.n), ik2(sym.n), ik3(sym.n);
complex iomega, iomega2, iomega3;

macro RRT(u) (exp(-paramTa/u#T)/u#T^1.5) // EOM
macro dRRT(u) ((paramTa/u#T - 1.5)/u#T^2.5*exp(-paramTa/u#T)) // EOM
macro d2RRT(u) (((paramTa/u#T - 5.0)*paramTa/u#T + 3.75)/u#T^3.5*exp(-paramTa/u#T)) // EOM
macro d3RRT(u) ((((paramTa/u#T - 10.5)*paramTa/u#T + 26.25)*paramTa/u#T - 13.125)/u#T^4.5*exp(-paramTa/u#T)) // EOM

real constfact = params["Dhc"]*params["Da"]*(1.0 + paramTs)/(params["Re"]*paramAs)*(min(1.0, params["phi"])/(params["phi"] + 2.0/paramvfO2))^1.5;
//constfact /= params["Dhc"]*min(1.0, params["phi"])/(params["phi"] + paramAFR);
macro Q() ( y*constfact*(RRY(ub)*RRT(ub)) ) // EOM
macro Qp(u) ( y*constfact*(dRRY(ub)*u#Y*RRT(ub) + RRY(ub)*dRRT(ub)*u#T) ) //EOM
macro Qpp(u, uu) ( y*constfact*(d2RRY(ub)*RRT(ub)*u#Y*uu#Y + dRRY(ub)*dRRT(ub)*(u#Y*uu#T + u#T*uu#Y) + RRY(ub)*d2RRT(ub)*u#T*uu#T) ) //EOM
macro Qppp(u, uu, uuu) ( y*constfact*(d3RRY(ub)*RRT(ub)*u#Y*uu#Y*uuu#Y + d2RRY(ub)*dRRT(ub)*(u#Y*uu#Y*uuu#T + (u#Y*uu#T + u#T*uu#Y)*uuu#Y) + dRRY(ub)*d2RRT(ub)*((u#Y*uu#T + u#T*uu#Y)*uuu#T + u#T*uu#T*uuu#Y) + RRY(ub)*d3RRT(ub)*u#T*uu#T*uuu#T) ) //EOM
  
  sym = 0;
  real R = int2d(Th)(Q);
  complex Rc = 0;
  for (int j = 0; j < Nh; j++){
    um[] = uh(:, j);
    um2[] = conj(um[]);
    ik.im = (1+j)*(sym0);
    ik2 = conj(ik);
    iomega = 1i*(1+j)*(omega);
    iomega2 = -iomega;
    Rc = int2d(Th)(Qpp(um, um2));
    R += real(Rc);
    for (int k = 0; k < Nh-j-1; k++){
      um2[] = uh(:, k);
      um3[] = conj(uh(:, 1+j+k));
      ik2.im = (1+k)*(sym0);
      ik3.im = -(2+j+k)*(sym0);
      iomega2 = 1i*(1+k)*(omega);
      iomega3 = -1i*(2+j+k)*(omega);
      Rc = int2d(Th)(Qppp(um, um2, um3));
      R += real(Rc);
    }
  }
cout.precision(17);
real Qbase = int2d(Th)(Q);
real Qmean = R;
cout << "Qbase = \n\t" << Qbase << endl;
cout << "Qmean = \n\t" << Qmean << endl;
```