%
%Author  : Ender Konukoglu
%Email   : ender.konukoglu@gmail.com
%Article : @inproceedings{konukoglu2007recursive,
%          title={A recursive anisotropic fast marching approach to reaction diffusion equation: Application to tumor growth modeling},
%          author={Konukoglu, Ender and Sermesant, Maxime and Clatz, Olivier and Peyrat, Jean-Marc and Delingette, Herve and Ayache, Nicholas},
%          booktitle={Information processing in medical imaging},
%          pages={687--699},
%          year={2007},
%          organization={Springer}
%					}
%Date    : June 21, 2013.
%
slab = 70;
boundary = zeros(slab,slab,slab);
object = zeros(slab,slab,slab); object(35,35,35) = 1;
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
figure,imagesc(Distance(:,:,35));
figure,imagesc(squeeze(Distance(:,35,:)));
figure,imagesc(squeeze(Distance(35,:,:)));
