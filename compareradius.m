ptsz = size(tree, 1);
[Iout,whatScale,Voutx,Vouty,Voutz]=FrangiFilter3D(I_original);
[M, N, Z] = size(I_original);
for i = 1 : ptsz
	input_x = round(tree(i, 3));
	input_y = round(tree(i, 4));
	input_z = round(tree(i, 5));
	whatScale(input_x, input_y, input_z)
	radius = getradiusfrangi(I_original, whatScale, input_x, input_y, input_z);
	tree(i, 6) = radius;
end
figure
hold on
for i = 1 : M
	for j = 1 : N
		for k = 1 : Z
			if  (whatScale(i, j, k) > 1)
				% plot3(i, j, k, '.', 'Color', [0.2 0.2 0.2]);
				switch whatScale(i, j, k)
				    case 2
				        plot3(i, j, k, '.', 'Color', [0.1 0.1 0.1]);
				    case 3
				        plot3(i, j, k, '.', 'Color', [0.2 0.2 0.2]);
				        % disp('Text');
				    case 4
				        plot3(i, j, k, '.', 'Color', [0.7 0.7 0.7]);
				end
				
			end
		end
	end
end
plot3(x, y, z, '.', 'Color', [0 0 1]);
hold off
