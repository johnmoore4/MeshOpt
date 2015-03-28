ff_len = 10;
max_lan = ff_len;
surf_len = 0.1;
dist = 10.0;
cuttof_dist = 1;

Geometry.Tolerance = 1.0e-12;

Mesh.CharacteristicLengthFactor=2.0;

Mesh.CharacteristicLengthMax = ff_len;

Mesh.SaveParametric = 1;

Merge "Sphere25CSymmetry.stp";


Physical Line(0) = {1:24};

Physical Surface(20) = {10};
Physical Surface(6) = {5:9};
Physical Surface(7) = {1:4};

Physical Volume(0) = {1};

Field[5] = Box;
Field[5].VOut = ff_len;
Field[5].VIn = surf_len;
Field[5].XMax = 1;
Field[5].XMin = -1;
Field[5].YMax = 1;
Field[5].YMin = -1;
Field[5].ZMax = 1;
Field[5].ZMin = -1;

Field[6] = Box;
Field[6].VOut = ff_len;
Field[6].VIn = 0.25;
Field[6].XMax = 10;
Field[6].XMin = -1.5;
Field[6].YMax = 2;
Field[6].YMin = -2;
Field[6].ZMax = 2;
Field[6].ZMin = -2;


Field[8] = Box;
Field[8].VOut = ff_len;
Field[8].VIn = surf_len;
Field[8].XMax = 0.75;
Field[8].XMin = -0.75;
Field[8].YMax = 0.75;
Field[8].YMin = -0.75;
Field[8].ZMax = 0.75;
Field[8].ZMin = -0.75;

Field[10] = Min;
Field[10].FieldsList = {5,6,8};

Background Field = 10;

