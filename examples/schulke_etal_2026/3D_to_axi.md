# 3D_to_axi.md
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
assert(mpisize==1);
string fileroot, fileext = parsefilename(filein, fileroot);
parsefilename(fileout, fileout); // trim extension from output file, if given
if(filein != "" && meshin == "") meshin = readmeshname(workdir + filein); // get mesh file
string meshroot, meshext = parsefilename(meshin, meshroot);
cout << meshroot << fileroot << endl;
// Load mesh, make FE basis
Th = readmeshN(workdir + meshin);
Thg = Th;
DmeshCreate(Th);
restu = restrict(XMh, XMhg, n2o);
XMh defu(ub);
if(fileext == "base") {
    if (mpirank == 0) cout << "  Loading '" + fileroot + ".base' on '" + meshin + "' from '" + workdir + "'." << endl;
    string filemesh;
    ifstream file(workdir + fileroot + ".base");
    file >> filemesh >> filemesh;
    if(paramnames[0] != ""){
    	for (int k = 0; k < paramnames.n; ++k)
        	file >> paramnames[k] >> params[paramnames[k]];}
    if(monitornames[0] != ""){
	  	for (int k = 0; k < monitornames.n; ++k)
			file >> monitornames[k] >> monitors[monitornames[k]];}
    if (filemesh == meshin){ // no interpolation needed
        fespace XMhg1(Thg, [P2, P2, P2, P2, P2, P1]);
        XMhg1 [ubg1, ubg1y, ubg1z, ubg1Y, ubg1T, ubg1p];
        file >> ubg1[];
		defu(ub) = [ubg1, ubg1y, ubg1Y, ubg1T, ubg1p];
    }
    else { // must interpolate
        if (mpirank == 0) cout << "\tMesh mismatch. Interpolating from '" << filemesh << "'." << endl;
        meshN Thg1 = readmeshN(workdir + filemesh);
        fespace XMhg1(Thg1, [P2, P2, P2, P2, P2, P1]);
        XMhg1 [ubg1, ubg1y, ubg1z, ubg1Y, ubg1T, ubg1p];
        file >> ubg1[];
		defu(ub) = [ubg1, ubg1y, ubg1Y, ubg1T, ubg1p];
    }
    setparams(paramnames, params);
}
savebase(fileout, "", meshin, true, true);
```