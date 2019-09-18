p1 = newp; Point(p1) = {-0.0206, 0, -0.0005};
p2 = newp; Point(p2) = {-0.0206, 0.0254, -0.0005};
p3 = newp; Point(p3) = {0, -0.0254, -0.0005};
p4 = newp; Point(p4) = {0, 0, -0.0005};
p5 = newp; Point(p5) = {0, 0.0254, -0.0005};
p6 = newp; Point(p6) = {0.206, -0.0254, -0.0005};
p7 = newp; Point(p7) = {0.206, 0.007, -0.0005};
p8 = newp; Point(p8) = {0.206, 0.0254, -0.0005};
p9 = newp; Point(p9) = {0.290, -0.0166, -0.0005};
p10 = newp; Point(p10) = {0.290, 0, -0.0005};
p11 = newp; Point(p11) = {0.290, 0.0166, -0.0005};
Line(1) = {5, 2};
Line(2) = {2, 1};
Line(3) = {1, 4};
Line(4) = {4, 5};
Line(5) = {4, 3};
Line(6) = {3, 6};
Line(7) = {6, 7};
Line(8) = {7, 4};
Line(9) = {5, 8};
Line(10) = {8, 7};
Line(11) = {6, 9};
Line(12) = {9, 10};
Line(13) = {10, 7};
Line(14) = {10, 11};
Line(15) = {11, 8};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Curve Loop(2) = {9, 10, 8, 4};
Plane Surface(2) = {2};
Curve Loop(3) = {5, 6, 7, 8};
Plane Surface(3) = {3};
Curve Loop(4) = {15, 10, -13, 14};
Plane Surface(4) = {4};
Curve Loop(5) = {7, -13, -12, -11};
Plane Surface(5) = {5};

prog = 1.1;
middle1 = 200;
middle2 = 60;
yn = 40;
bumpRatio = 0.05;


Transfinite Curve{-1,3} = 20 Using Progression 0.95;
Transfinite Curve{4,-2} = yn Using Bump bumpRatio; //Progression prog;
Transfinite Surface{1} = {1,4,5,2};

Transfinite Curve{-8,9} = middle1 Using Progression 1.005;
Transfinite Curve{-10,4} = yn Using Bump bumpRatio; //Progression prog;
Transfinite Surface{2} = {4,7,8,5};

Transfinite Curve{-8,6} = middle1 Using Progression 1.005;
Transfinite Curve{5,-7} = yn Using Bump bumpRatio; //Progression prog;
Transfinite Surface{3} = {6,7,4,3};

Transfinite Curve{-13,-15} = middle2 Using Progression 1;
Transfinite Curve{-10,14} = yn Using Bump bumpRatio; //Progression prog;
Transfinite Surface{4} = {10,11,8,7 };

Transfinite Curve{13,11} = middle2 Using Progression 1;
Transfinite Curve{-7,-12} = yn Using Bump bumpRatio; //Progression prog;
Transfinite Surface{5} = {9,10,7,6};


Recombine Surface{1};
Recombine Surface{2};
Recombine Surface{3};
Recombine Surface{4};
Recombine Surface{5};//+
Extrude {0, 0, 0.001} {
  Surface{2}; Surface{1}; Surface{3}; Surface{4}; Surface{5}; Layers{1}; Recombine;
}
//+
Physical Volume("body") = {2, 1, 3, 4, 5};
//+
Physical Surface("inlet") = {50};
//+
Physical Surface("outlet") = {102, 120};
//+
Physical Surface("upperWall") = {46, 24, 90};
Physical Surface("lowerWall") = {68, 72, 124};
