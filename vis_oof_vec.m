% This script tries to visualise eigen vector direction of optimal oriented flux 
load('mat\diffussion.mat');
load('mat\bI.mat');
nx = size(I, 1);
ny = size(I, 2);
nz = size(I, 3);
figure
hold on
fprintf('highlight the direction of first eigenvector\n');
scale_coefficient = 10;
for i = 1 : nx
    for j = 1 : ny
        for k = 1 : nz
			if I(i,j,k) == 1
                vecI = squeeze(T(i,j,k,:));
                hessianmat = hessianvaluetomat(vecI);
                [V, D] = eig(hessianmat);
                tensormat = scale_coefficient * D(1,1) * V(:,1) * V(:,1)' + D(2,2) * V(:,2) * V(:,2)' + D(3,3) * V(:,3) * V(:,3)';
                newvecI = [tensormat(1,1); tensormat(1,2); tensormat(1,3); tensormat(2,2); tensormat(2,3); tensormat(3,3);]; 
                newT(i,j,k,:) = newvecI; 
				% quiver3(i, j, k, V(1,1), V(2,1), V(3,1));
            else                
                newT(i,j,k,:) = squeeze(T(i,j,k,:));
            end
        end
    end
end
hold off
plot_hessian_whole(newT, I)
% A * V - V * D;
% T = D(1,1) * V(:,1) * V(:,1)' + D(2,2) * V(:,2) * V(:,2)' + D(3,3) * V(:,3) * V(:,3)';