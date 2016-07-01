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
slab = 512;
boundary = zeros(slab,slab);
object = zeros(slab,slab); object(slab/2,slab/2) = 1;
volDim = [1,1,1];
F = ones(slab,slab) * 2;
angle = pi/7;
V = [sin(angle), cos(angle); -cos(angle), sin(angle)];
L = [0.3, 0;0, 1.2;];
tensor = V * L * V';
T = zeros(slab,slab,3);
T(:,:,1) = ones(slab,slab) * tensor(1,1); 
T(:,:,2) = ones(slab,slab) * tensor(1,2);
T(:,:,3) = ones(slab,slab) * tensor(2,2);
tic;
Distance = mxAnisoDistanceTransform(object, T, boundary, F, volDim);
toc
figure,imagesc(squeeze(Distance));
