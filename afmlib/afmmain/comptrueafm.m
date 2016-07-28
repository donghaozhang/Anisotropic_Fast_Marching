slab = 20;
boundary = zeros(slab,slab,slab);
object = zeros(slab,slab,slab); object(round(slab/2),round(slab/2),round(slab/2)) = 1;
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
tic;
Distance = mxAnisoDistanceTransform(object, T, boundary, F, volDim);
toc
Distance = permute(Distance, [3 2 1])
figure,imagesc(Distance(:,:,round(slab/2)));
figure,imagesc(squeeze(Distance(:,round(slab/2),:)));
figure,imagesc(squeeze(Distance(round(slab/2),:,:)));