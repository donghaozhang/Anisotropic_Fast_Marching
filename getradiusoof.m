function radius = getradiusoof(inputMatrix, input_x, input_y, input_z, eigoofmat)
	[M, N, Z] = size(inputMatrix);
	input_x = constrain(input_x, 1, M);
	input_y = constrain(input_y, 1, N);
	input_z = constrain(input_z, 1, Z);
    input_x = round(input_x);
    input_y = round(input_y);
    input_z = round(input_z);
	radius = eigoofmat(input_x, input_y, input_z, 1) + eigoofmat(input_x, input_y, input_z, 2) + eigoofmat(input_x, input_y, input_z, 3);
	% radius = eigoofmat(input_x, input_y, input_z, 1);
	radius = ceil(radius) + 3;
	if radius < 0
		radius = 1;
	end
end