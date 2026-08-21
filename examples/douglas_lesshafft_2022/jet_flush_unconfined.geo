//
// jet_flush_unconfined.geo
// Chris Douglas
// christopher.douglas@duke.edu
//
// This file can be used with FreeFEM to create a mesh for the flush-mounted, 
// radially-unconfined jet as in [Douglas & Lesshafft. JFM, (2022)].
//

n0 = 32;
n1 = 8;
n2 = 1;
l = 4.0;
deltaw = 1.0e-5;
xw = 1.0;
Rw = 8.0;
Rinf = 40.0;
Xinf = 100.0;

// Points
//          4---3--..,_
//          |   |      `-.
//          |   |         `.
//          5---6           `
//  8-----------7            \
//  |                         |
//  1-----------0-------------2
Point (0)  = {0, 0, 0, 1/n1}; // origin
Point (1)  = {-l, 0, 0, 1/n1};
Point (2)  = {Rinf, 0, 0, 1/n2};
Point (3)  = {0, Rinf, 0, 1/n2};
//Point (4)  = {-xw, Rinf, 0, 1/n2};
//Point (5)  = {-xw, 0.5+deltaw, 0, 1/n1};
Point (6)  = {0, 0.5+deltaw, 0, 1/n0};
Point (7)  = {0, 0.5, 0, 1/n0};
Point (8)  = {-l, 0.5, 0, 1/n1};

// Lines
Line (1) = {8,1};
Line (2) = {1,2};
Circle (3) = {2,0,3};
//Line (4) = {3,4};
//Line (5) = {4,5};
//Line (6) = {5,6};
Line (7) = {6,7};
Line (8) = {7,8};
Line (9) = {3,6};

// Labels
Physical Line ("AXIS") = {2};
Physical Line ("OUTFLOW") = {3};
Physical Line ("INFLOW") = {1};
Physical Line ("WALL") = {9};
Physical Line ("PIPE") = {8};

// Surfaces
Line Loop(1) = {1:3,7:9};
Plane Surface(1) = {1};
Physical Surface("DOMAIN") = {1};
