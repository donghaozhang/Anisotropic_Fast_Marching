clear all
clc
close all
slab = 10;
boundary = zeros(slab,slab,slab);
object = zeros(slab,slab,slab); object(5,5,5) = 1;
volDim = [1,1,1];
F = ones(slab,slab,slab) * 2;
angle = pi/7;
V = [sin(angle), cos(angle), 0; -cos(angle), sin(angle), 0; 0, 0, 1];
L = [0.3, 0, 0;0, 10, 0;0, 0, 3];
tensor = V * L * V';
T = zeros(slab,slab,slab,6);
T(:,:,:,1) = ones(slab,slab,slab) * tensor(1,1); 
T(:,:,:,2) = ones(slab,slab,slab) * tensor(1,2);
T(:,:,:,3) = ones(slab,slab,slab) * tensor(1,3);
T(:,:,:,4) = ones(slab,slab,slab) * tensor(2,2);
T(:,:,:,5) = ones(slab,slab,slab) * tensor(2,3);
T(:,:,:,6) = ones(slab,slab,slab) * tensor(3,3);
afmSize = [size(T, 1) size(T, 2) size(T, 3);];
T = permute(T, [4 3 2 1]);
F = permute(F, [3 2 1]);
boundary = permute(boundary, [3 2 1]);
FO = permute(object, [3 2 1]);
dx = 1; dy = 1; dz = 1;
timeLimit = 10000;

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
[Tval, Ttag] = afmAnisotropicFastMarching(Tval, Ttag, F, boundary, dx, dy, dz, afmSize, T, timeLimit);
Tval(Tval==100) = -1;
figure,imagesc(Tval(:,:,5));
figure,imagesc(squeeze(Tval(:,5,:)));
figure,imagesc(squeeze(Tval(5,:,:)));