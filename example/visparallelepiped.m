figure
axis equal;
hold on
X = 0;
Y = 0;
Z = 0;
U = 1;
V = 1;
W = 1;
quiver3(X,Y,Z,U,V,W,0.5)
U = 1;
V = -1;
W = 1;
quiver3(X,Y,Z,U,V,W,0.5)
U = 1;
V = 1;
W = -1;
quiver3(X,Y,Z,U,V,W,0.5)
vecA = [1, 1, 1];
vecB = [1, -1, 1];
vecC = [1, 1, -1];

BcrossC = cross(vecB, vecC);
AdotBC = dot(vecA, BcrossC);