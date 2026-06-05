//
// vortex.geo
// Chris Douglas
// christopher.douglas@duke.edu
//
// This file can be used with Gmsh to create a mesh for the Grabowski-Berger
// vortex problem.
//

n0 = 10;
n1 = 3;
n2 = 1;
xmax = 120.0;
rmax = 60.0;

// Points
//  3-------------------------2
//  |                         |
//  |                         |
//  |                         |
//  |                         |
//  0-------------------------1
Point (0)  = {0, 0, 0, 1/n0}; // origin
Point (1)  = {xmax, 0, 0, 1/n2};
Point (2)  = {xmax, rmax, 0, 1/n2};
Point (3)  = {0, rmax, 0, 1/n1};


// Lines
Line (1) = {0,1};
Line (2) = {1,2};
Line (3) = {2,3};
Line (4) = {3,0};

// Labels
Physical Line ("AXIS") = {1};
Physical Line ("OUTFLOW") = {2};
Physical Line ("LATERAL") = {3};
Physical Line ("INFLOW") = {4};

// Surfaces
Line Loop(1) = {1:4};
Plane Surface(1) = {1};
Physical Surface("DOMAIN") = {1};
