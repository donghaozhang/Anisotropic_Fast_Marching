% The correct spelling of diffusion should be diffusion rather than diffussion.
load('../../mat/diffussion.mat');
x = 100; 
y = 32;
z = 14;
D = permute(T, [4 3 2 1]);
F = 40000000;
yon = [1, -4, 9]; 
yon = yon';
detDmat = afmgroup_velocity(x, y, z, D, yon, F)