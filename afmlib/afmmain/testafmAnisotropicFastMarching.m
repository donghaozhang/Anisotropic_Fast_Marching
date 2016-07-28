	% inputs:
	% BI     : boundary volume
	% T      : tensor volume (D)
	% FO     : first object
	% F      : speed image
	% volDim : voxel dimensions
	% outputs: 
	% OI : distance transform
clc;
clear all;
close all;
load('../../mat/diffussion.mat');
afmSize = [size(T, 1) size(T, 2) size(T, 3);];
F = ones(afmSize(3), afmSize(2), afmSize(1));
Boundary = zeros(afmSize(3), afmSize(2), afmSize(1));
dx = 1;
dy = 1;
dz = 1;
timeLimit = 10000;

% Constructing Volumes 
fprintf('Constructing volumes...\n');
D =  permute(T, [4 3 2 1]);
FO = zeros(afmSize(3), afmSize(2), afmSize(1));
FO(20, 55, 122) = 1;  
for(z = 1 : afmSize(3))
	% fprintf('after constructing volumes Current z value: %d\n', z);
	for(y = 1 : afmSize(2))
		for(x = 1 : afmSize(1))
			if (FO(z, y, x) == 1)
				Tval(z,y,x) = 0;
				Ttag(z,y,x) = 200;
			else
				Tval(z,y,x) = 100;
				Ttag(z,y,x) = 100;
			end
        end
    end
end
[Tval, Ttag] = afmAnisotropicFastMarching(Tval, Ttag, F, Boundary, dx, dy, dz, afmSize, D, timeLimit);
Tval(Tval==100) = -1;

	% /* Assigning the output */
	% for (int z = 0; z < dimV[2]; z++)
	% 	for (int y = 0; y < dimV[1]; y++)
	% 		for (int x = 0; x < dimV[0]; x++)
	% 		{
	% 			index = x + y * dimV[0] + z * dimV[0] * dimV[1];
	% 			if (Ttag[z][y][x] == 200 && (unsigned char)FO[index] == 0)
	% 				OI[index] = (double)Tval[z][y][x];
	% 			else
	% 				OI[index] = -1.0f;
	% 		}