% in this script I am going to plot the direction of vessel direction using using quiver3
load('zebraI.mat');
load('mat\bI.mat');
[Iout,whatScale,Voutx,Vouty,Voutz]=FrangiFilter3D(I);
nx = size(I, 1);
ny = size(I, 2);
nz = size(I, 3);
figure
hold on 
for i = 1 : nx
    for j = 1 : ny
        for k = 1 : nz
            if I(i,j,k) == 1
				quiver3(i,j,k,Voutx(i,j,k),Vouty(i,j,k),Voutz(i,j,k));
            end
        end
    end
end
hold off
