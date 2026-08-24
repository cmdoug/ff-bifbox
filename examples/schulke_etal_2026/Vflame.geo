//
// Vflame.geo
// Chris Douglas
// christopher.douglas@duke.edu
//
// This file can be used with Gmsh to create a mesh for the annular V-flame configuration.
// It is based on a similar file from `ff-bifbox/examples/wang_etal_2024/Vflame.md`.
//

n0 = 125;
n1 = 31.25;
n2 = 25;
n3 = 12.5;
n4 = 5;

xmin = -3;
xmax = 20;
rmax = 7.5;
D = 1.;
Dref = 1.1;
Dcb = 0.3/Dref;
xcb = 0.2/Dref;

right_end = 3.8/Dref;
up_in = 0.8/Dref;
up_cb = 1.0/Dref;
up_middle = 2.5/Dref;
up_end = 3.5/Dref;
left_in = -0.5/Dref;
x1 = 6.0/Dref;
r1 = 4.0/Dref;

// Points
//             22----------------------------------------21
//              |                                        |
//             19-------------------18                   |
//              |                   |                    |
//              9-----------11      |                    |
//              |           |       |                    |
//              | ----------12      |                    |
//              |/          |       |                    |
//              8           |       |                    |
//              |           |       |                    |
//  4-----------7           |       |                    |
//  |                       |       |                    |
//  |           ------------15      |                    |
//  |          /            |       |                    |
//  |         /      -------16      |                    |
//  3--------14---2 /       |       |                    |
//                |/        |       |                    |
//              o 1---------13------17-------------------20
Point(1) = {xcb, 0, 0, 1/n0};
Point(2) = {xcb, 0.5*Dcb, 0, 1/n0};
Point(3) = {xmin, 0.5*Dcb, 0, 1/n2};
Point(4) = {xmin, 0.5*D, 0, 1/n2};
Point(7) = {0, 0.5*D, 0, 1/n2};
Point(8) = {0, up_in, 0, 1/n2};
Point(9) = {0, up_end, 0, 1/n2};
Point(11) = {right_end, up_end, 0, 1/n2};
Point(12) = {right_end, up_middle, 0, 1/n1};
Point(13) = {right_end, 0, 0, 1/n2};
Point(14) = {left_in, 0.5*Dcb, 0, 1/n0};
Point(15) = {right_end, up_cb+up_in, 0, 1/n0};
Point(16) = {right_end, up_cb, 0, 1/n0};
Point(17) = {x1, 0, 0, 1/n3};
Point(18) = {x1, r1, 0, 1/n3};
Point(19) = {0, r1, 0, 1/n3};
Point(20) = {xmax, 0, 0, 1/n4};
Point(21) = {xmax, rmax, 0, 1/n4};
Point(22) = {0, rmax, 0, 1/n4};

// Lines
Line(1) = {3, 14};
Line(2) = {14, 2};
Line(3) = {2, 1};
Line(4) = {1, 13};
Line(5) = {13, 16};
Line(6) = {16, 15};
Line(7) = {15, 12};
Line(8) = {12, 11};
Line(9) = {11, 9};
Line(10) = {9, 8};
Line(11) = {8, 7};
Line(12) = {7, 4};
Line(13) = {4, 3};
Line(15) = {14, 15};
Line(16) = {1, 16};
Line(17) = {12, 8};
Line(18) = {19, 9};
Line(19) = {13, 17};
Line(20) = {17, 18};
Line(21) = {18, 19};
Line(22) = {19, 22};
Line(23) = {22, 21};
Line(24) = {21, 20};
Line(25) = {20, 17};

// Labels
Physical Line("axis") = {4, 19, 25};
Physical Line("outlet") = {24};
Physical Line("lateral") = {23};
Physical Line("planewall") = {10, 11, 18, 22};
Physical Line("annuluswall") = {1, 2, 12};
Physical Line("inlet") = {13};
Physical Line("centerbodywall") = {3};

// Surfaces
Line Loop(1) = {-23, -24, -25, 19, 4, 3, 2, 1, 13, 12, 11, 10, 18, -22};
Plane Surface(1) = {1};
Line {5:9,15:17,20:21} In Surface{1};
Physical Surface("DOMAIN") = {1};