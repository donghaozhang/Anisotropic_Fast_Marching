function radius = getradiusfrangi(inputMatrix, responsematrix, input_x, input_y, input_z)
	%Calculate the radius of M * N * Z binary matrix at specific location 
	%The specific location is defined by corresponding  input_x, input_y, input_z
	% responsematrix uses the radius which can generate maximum response of frangi
	[M, N, Z] = size(inputMatrix);
	input_x = round(input_x);
	input_y = round(input_y);
	input_z = round(input_z);
	input_x = constrain(input_x, 1, M);
	input_y = constrain(input_y, 1, N);
	input_z = constrain(input_z, 1, Z);
	radius = responsematrix(input_x, input_y, input_z);
	radius = radius + 1;
end