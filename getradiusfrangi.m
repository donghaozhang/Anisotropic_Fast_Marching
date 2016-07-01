function radius = getradiusfrangi(inputMatrix, responsematrix, input_x, input_y, input_z)
	%Calculate the radius of M * N * Z binary matrix at specific location 
	%The specific location is defined by corresponding  input_x, input_y, input_z
	% responsematrix uses the radius which can generate maximum response of frangi
	[M, N, Z] = size(inputMatrix);
	input_x = constrain(1, M);
	input_y = constrain(1, N);
	input_z = constrain(1, Z);
	radius = responsematrix(input_x, input_y, input_z);
end