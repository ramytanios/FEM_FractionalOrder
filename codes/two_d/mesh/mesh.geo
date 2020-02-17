//+
SetFactory("Built-in");
//+
SetFactory("Built-in");
//+
SetFactory("Built-in");
//+
SetFactory("Built-in");
//+
SetFactory("Built-in");
//+
SetFactory("Built-in");
//+
Point(1) = {1, 1, 0, 1.0};
//+
Point(2) = {1, 8, 0, 1.0};
//+
Point(3) = {8, 8, 0, 1.0};
//+
Point(4) = {8, 1, 0, 1.0};
//+
Point(5) = {8, 1, 0, 1.0};
//+
Line(1) = {2, 3};
//+
Line(2) = {4, 3};
//+
Line(3) = {2, 1};
//+
Line(4) = {1, 4};
//+
Physical Curve("e11") = {1};
//+
Physical Curve("e22") = {2};
//+
Physical Curve("e33") = {4};
//+
Physical Curve("e44") = {3};
//+
Curve Loop(1) = {1, -2, -4, -3};
//+
Plane Surface(1) = {1};
//+
Physical Surface("surf") = {1};
//+
Transfinite Surface {1};
//+
Transfinite Surface {1};
