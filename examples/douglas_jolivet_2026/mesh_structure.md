# mesh_structure.md
Author: Chris Douglas ([@cmdoug](https://github.com/cmdoug)) [christopher.douglas@duke.edu](mailto:christopher.douglas@duke.edu)

This file is used to create initial meshes for the example of the buckling 3-D structure from [C. Douglas, P. Jolivet. (2026)].

```freefem
assert(mpisize == 1); // Must be run with 1 processor
include "settings.idp" // so that boundary labels match the convention in the settings file
real R = 2.540; //[m]
real t = 0.00635; //[m]
real beta = 0.1; //[rad]
int n0 = 2*ceil(getARGV("-n0", 40)/2); //must be even for discrete symmetry and point loading
int n1 = 2*ceil(getARGV("-n1", 2)/2); //must be even for pinned BC
int flag = getARGV("-flag", 1);
string meshroot, meshout = getARGV("-mo", "shell.mesh"); // mesh filename

func string parsefilename(string & filename, string & fileroot){
  string fileext;
  if(filename.rfind(".") > 0){ // filename includes extension
    fileext = filename(filename.rfind(".")+1:filename.length-1); // get file extension
    fileroot = filename(0:filename.rfind(".")-1); // get file root
  }
  else {
    fileroot = filename;
    fileext = "";
  }
  return fileext;
}

if(meshout.rfind(".mesh") < 0) meshout = meshout + ".mesh"; // add extension if not provided
string meshext = parsefilename(meshout, meshroot);

string[int] symlist = ["S", "Hx", "Hy", "Rz", "A"];

int[int] rup = [0, BCtop], rdown = [0, BCfree], rmid(8);
mesh Th = square(n0/2, n0/2, [beta*(x - 1.), beta*(y - 1.)], flags = flag, region = 0);
int[int] regnum = [0, 0];
mesh3 Th3 = buildlayers(Th, n1, zbound = [-t/2., t/2.], region = regnum, labelup = rup, labeldown = rdown);
for (int ii = 0; ii < 5; ii++) {
  string meshsym = symlist[ii];
  if (meshsym == "Rz") regnum = [0, BCpsym];
  else regnum = [0, 0];
  mesh3 Thg = movemesh(Th3, [R*x, (R+z)*sin(y), (R+z)*cos(y) - R*cos(beta)], region = regnum);
  if (meshsym == "S") rmid = [1, BCpin, 2, BCxsym, 3, BCysym, 4, BCfree];
  else if (meshsym == "Hx") {
    mesh3 Thy = movemesh(Thg, [x, -y, z], orientation = -1);
    Thg = Thg + Thy;
    rmid = [1, BCpin, 2, BCxsym, 3, BCfree, 4, BCfree];
  }
  else if (meshsym == "Hy" || meshsym == "A" || meshsym == "Rz") {
    mesh3 Thx = movemesh(Thg, [-x, y, z], orientation = -1);
    Thg = Thg + Thx;
    if (meshsym == "A" || meshsym == "Rz"){
      Thx = movemesh(Thg, [x, -y, z], orientation = -1);
      Thg = Thg + Thx;
      rmid = [1, BCpin, 2, BCfree, 3, BCfree, 4, BCfree];
    }
    else rmid = [1, BCpin, 2, BCfree, 3, BCysym, 4, BCfree];
  }
  else exit(1);
  Thg = change(Thg, label = rmid, rmInternalFaces = true);
  int[int] meshlabels = labels(Thg);
  cout << "\tMesh: " << Thg.nv << " vertices, " << Thg.nt << " elements, " << Thg.nbe << " boundary elements, " << meshlabels.n << " labeled boundaries." << endl;
  cout << "  Saving mesh '" + meshroot + "_" + meshsym + "." + meshext + "'." << endl;
  savemesh(Thg, meshroot + "_" + meshsym + "." + meshext);
}
```