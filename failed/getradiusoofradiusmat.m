function radius = getradiusoofradiusmat(inputMatrix, input_x, input_y, input_z, eigoofmat)
	[M, N, Z] = size(inputMatrix);
	input_x = constrain(input_x, 1, M);
	input_y = constrain(input_y, 1, N);
	input_z = constrain(input_z, 1, Z);
	radius = eigoofmat(input_x, input_y, input_z);
end