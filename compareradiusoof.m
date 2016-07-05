load('mat\eigvoof.mat');
% load('mat\radius_output.mat');
ptsz = size(tree, 1);
for i = 1 : ptsz
	input_x = round(tree(i, 3));
	input_y = round(tree(i, 4));
	input_z = round(tree(i, 5));
	radius = getradiusoof(I_original, input_x, input_y, input_z, eigvoof);
	% radius = getradiusoofradiusmat(I_original, input_x, input_y, input_z, radius_output);
	% tree(i, 6) = (radius + tree(i, 6)) / 2;
	% tree(i, 6) = radius / 3;
	tree(i, 6) = radius;
end