
fflc = 10.0;
surf_len = 0.1;
lete_factor = 0.25;

Geometry.Tolerance = 1.0e-12;
Mesh.CharacteristicLengthMax = fflc;
Mesh.CharacteristicLengthExtendFromBoundary = 1;
Mesh.SaveParametric = 1;
Merge "rae2822Square.stp";

// Set physical tags for boundary conditions
Physical Line(6) = {1:4};
Physical Line(7) = {5,6};


Physical Surface(0) = {1};


Field[1] = Attractor;
Field[1].NodesList = {5,6};


Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = 0.01;
Field[2].LcMax = 0.05;
Field[2].DistMin = 0.0;
Field[2].DistMax = 0.15;
Field[2].StopAtDistMax = 1;


Field[5] = Box;
Field[5].VOut = fflc;
Field[5].VIn = surf_len;
Field[5].XMax = 5;
Field[5].XMin = -1;
Field[5].YMax = 1.0;
Field[5].YMin = -1.0;
Field[5].ZMax = 1;
Field[5].ZMin = -1;

Field[6] = Box;
Field[6].VOut = fflc;
Field[6].VIn = 5*surf_len;
Field[6].XMax = 25.0;
Field[6].XMin = -5;
Field[6].YMax = 3.0;
Field[6].YMin = -3.0;
Field[6].ZMax = 1;
Field[6].ZMin = -1;


Field[10] = Min;
Field[10].FieldsList = {2,5,6};

Background Field = 10;



