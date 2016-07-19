load('../../mat/diffussion.mat');
x = 100; 
y = 32;
z = 14;
dx = 1;
dy = 1;
dz = 1;
tetNo = 8;
[trX, trY, trZ] =  afmSetTetrahedra();
D =  permute(T, [4 3 2 1]);
afmSize = [size(T, 1) size(T, 2) size(T, 3);];
Boundary = zeros(afmSize);
Boundary = permute(Boundary, [3 2 1]);
Tvalue = zeros(afmSize);
Tvalue = permute(Tvalue, [3 2 1]);
Ttag = zeros(afmSize);
Ttag = permute(Ttag, [3 2 1]);
[Tvalue, Ttag, value_iter] = afmUpdateNeighborhoodTrial(Tvalue, Ttag, Boundary, dx, dy, dz, afmSize, D, x, y, z, trX, trY, trZ, tetNo)