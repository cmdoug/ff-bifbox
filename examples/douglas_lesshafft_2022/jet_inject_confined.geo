//
// jet_inject_confined.geo
// Chris Douglas
// christopher.douglas@duke.edu
//
// This file can be used with FreeFEM to create a mesh for the axially-injected, 
// radially-confined jet as in [Douglas & Lesshafft. JFM, (2022)].
//

n0 = 32;
n1 = 8;
n2 = 1;
l = 4.0;
deltaw = 1.0e-5;
Xw = 1.0;
Rw = 8.0;
Rinf = 40.0;
Xinf = 100.0;

// Points
//          4---3-------------21
//          |   |             |
//          |   |             |
//          5---6-------------20
//  8-----------7             |
//  |                         |
//  1-----------0-------------2
Point (0)  = {0, 0, 0, 1/n1}; // origin
Point (1)  = {-l, 0, 0, 1/n1};
Point (2)  = {Xinf, 0, 0, 1/n2};
Point (20)  = {Xinf, 0.5+deltaw, 0, 1/n2};
Point (21)  = {Xinf, Rw, 0, 1/n2};
Point (3)  = {0, Rw, 0, 1/n2};
Point (4)  = {-Xw, Rw, 0, 1/n2};
Point (5)  = {-Xw, 0.5+deltaw, 0, 1/n1};
Point (6)  = {0, 0.5+deltaw, 0, 1/n0};
Point (7)  = {0, 0.5, 0, 1/n0};
Point (8)  = {-l, 0.5, 0, 1/n1};

// Lines
Line (1) = {8,1};
Line (2) = {1,2};
Line (20) = {2,20};
Line (21) = {20,21};
Line (3) = {21,3};
Line (4) = {3,4};
Line (5) = {4,5};
Line (6) = {5,6};
Line (7) = {6,7};
Line (8) = {7,8};
Line (9) = {3,6};
Line (10) = {6,20};

// Labels
Physical Line ("AXIS") = {2};
Physical Line ("OUTFLOW") = {20:21};
Physical Line ("INFLOW") = {1};
Physical Line ("WALL") = {3:7};
Physical Line ("PIPE") = {8};
Physical Line ("INTERNAL") = {9,10};

// Surfaces
Line Loop(1) = {1:2,20:21,3:8};
Plane Surface(1) = {1};
Line{9,10} In Surface{1};
Physical Surface("DOMAIN") = {1};
