
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
L = [1000 0;0 0.01];
t = [0:0.05:5000];
slab = 128;
d11_ring2 = ones(slab)*0.0001;
d12_ring2 = zeros(slab);
d22_ring2 = ones(slab)*0.0001;
for(k=1:length(t))

    x = [ -t(k)^1*sin(t(k)) + 1*t(k)^0.5*cos(t(k)), t(k)^1*cos(t(k)) + 1*t(k)^0.5*sin(t(k)) ]';    
    x = x/(norm(x) + 1e-5);
    y = [ t(k)^1*cos(t(k)) + 1*t(k)^0.5*sin(t(k)), t(k)^1*sin(t(k)) - 1*t(k)^0.5*cos(t(k)) ]';    
    y = y/(norm(y) + 1e-5);
    V = [x,y];
    D_ring2 = V*L*V';
    if(round(t(k)^1*cos(t(k))) + slab/2 > 0 && round(t(k)^1*cos(t(k))) + slab/2 < slab+1 ...
        && round(t(k)^1*sin(t(k))) + slab/2 > 0 && round(t(k)^1*sin(t(k))) + slab/2 < slab+1)
        d11_ring2( round(t(k)^1*cos(t(k))) + slab/2, round(t(k)^1*sin(t(k))) + slab/2 ) = D_ring2(1,1);
        d12_ring2( round(t(k)^1*cos(t(k))) + slab/2, round(t(k)^1*sin(t(k))) + slab/2 ) = D_ring2(1,2);
        d22_ring2( round(t(k)^1*cos(t(k))) + slab/2, round(t(k)^1*sin(t(k))) + slab/2 ) = D_ring2(2,2);
    end;
end;

se = gaussian_structural_element_2D([1 0;0 1],[2 2]);
se = se/sum(se(:));
d11_ring2 = convn(d11_ring2,se,'same');
d12_ring2 = convn(d12_ring2,se,'same');
d22_ring2 = convn(d22_ring2,se,'same');
boundary = zeros(slab,slab);
object = zeros(slab,slab); object(slab/2,slab/4) = 1;
volDim = [1,1,1];
F = ones(slab,slab);
T = zeros(slab,slab,3);
T(:,:,1) = d11_ring2; 
T(:,:,2) = d12_ring2;
T(:,:,3) = d22_ring2;
tic;
Distance = mxAnisoDistanceTransform(object, T, boundary, F, volDim);
toc
figure,imagesc(squeeze(Distance));
