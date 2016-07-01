ptsz = size(tree, 1);
[Iout,whatScale,Voutx,Vouty,Voutz]=FrangiFilter3D(I_original);
for i = 1 : ptsz
	input_x = round(tree(i, 3));
	input_y = round(tree(i, 4));
	input_z = round(tree(i, 5));
	radius = getradiusfrangi(I_original, whatScale, input_x, input_y, input_z);
	tree(i, 6) = radius;
end