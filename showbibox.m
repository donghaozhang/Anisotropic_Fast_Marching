function showbibox(B)
	[x y z] = ind2sub(size(B), find(B));
	plot3(y, x, z, '.', 'Color', [0.7 0.7 1]);  view(3); axis equal;
end