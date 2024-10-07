//Inputs
//https://www.youtube.com/watch?v=E8nGNZk71fw&ab_channel=OpenFOAMTutorials


boxdim = 1.;
//gridsize sert a rien
gridsize = boxdim/4;
numTransfinitSize = 3;

//C'est ce parametre qui controle la basesize des elements.
Mesh.CharacteristicLengthMax = 0.4;
Mesh.CharacteristicLengthMin = 0.4;
//Ne pas oublier de "Smooth 2D" le resultat
//Et on peut raffiner avec "Refine by splitting"


Point(1) = {0,0,0,gridsize};
Point(2) = {boxdim/2,0,0,gridsize};
Point(3) = {boxdim/2,boxdim,0,gridsize};
Point(4) = {0,boxdim,0,gridsize};

Point(5) = {boxdim,0,0,gridsize/10};
Point(6) = {boxdim,boxdim,0,gridsize/10};

Line(7) = {1,2};
Line(8) = {2,3};
Line(9) = {3,4};
Line(10) = {4,1};

Line(11) = {2,5};
Line(12) = {5,6};
Line(13) = {6,3};

Line Loop(14) = {7,8,9,10};
Line Loop(15) = {11,12,13,-8};

Plane Surface(16) = 14;
Plane Surface(17) = 15;

//Make one square structured.
Transfinite Line{7,8,9,10} = numTransfinitSize;
Transfinite Surface{16};
Recombine Surface{16};

Physical Line("Boundary 1") = {7,8,9,10};
Physical Line("Boundary 2") = {11,12,13,8};
Physical Surface("Surface Rect") = {16};
Physical Surface("Surface Tri") = {17};


//Recombine Surface{17};

Mesh.ElementOrder = 1;
Mesh.Algorithm = 5;

//Une fois le maillage generee, faire 
//"Smooth 2D"
//Et raffiner pour avoir une convergence en maillage