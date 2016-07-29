clear all
clc
close all
load('../../mat/prepareforafm.mat');
dx = 1; dy = 1; dz = 1;
timeLimit = 10000;
afmSize = [size(T, 1) size(T, 2) size(T, 3);];
T = permute(T, [4 3 2 1]);
FO = permute(object, [3 2 1]);
SpeedImage = permute(SpeedImage, [3 2 1]);
boundary_value = permute(boundary_value, [3 2 1]);
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
tic
[Tval, Ttag] = afmAnisotropicFastMarching(Tval, Ttag, SpeedImage, boundary_value, dx, dy, dz, afmSize, T, timeLimit);
Tval(Tval==100) = -1;