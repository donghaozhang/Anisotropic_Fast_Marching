function T_map = afm(I, threshold, foreground_speed_coefficient)
    bI = I > threshold;
    fprintf('Prepare for the speed image.\n');
    disp('Distance transform');
    notbI = not(I>threshold);
    bdist = bwdist(notbI, 'Quasi-Euclidean');
    bdist = bdist .* double(bI);
    bdist = double(bdist);
    [SourcePoint, maxD] = maxDistancePoint(bdist, I, true);
    % Speical treatment for anisotropic fast marching
    % SpeedImage= (bdist/maxD).^4;
    
    SpeedImage= (bdist/maxD) * foreground_speed_coefficient;
    background_speed = 1;
    SpeedImage(SpeedImage==0) = background_speed;

    % Original 
    % SpeedImage=(bdist/maxD).^4;
    % clear bdist;
    % SpeedImage(SpeedImage==0) = 1e-10;

    sigma_value = 0.43;
    % SpeedImage = double(bI) * 2;

    % figure(1), imagesc(squeeze(max(SpeedImage,[],3))), title('speed');

    % imagesc(max(DX,[],3))
    [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(double(I), sigma_value);
    szI = size(I);
    T = zeros(szI(1),szI(2),szI(3),6);

    eps = 1e-5;
    T(:,:,:,1) = Dxx;
    T(:,:,:,2) = Dxy;
    T(:,:,:,3) = Dxz;
    T(:,:,:,4) = Dyy;
    T(:,:,:,5) = Dyz;
    T(:,:,:,6) = Dzz;
    % why I do the following code is make sure that Dxx = 1; Dyy = 1; Dzz = 1;
    % Dxy == Dyx = 0 Dxz == Dzx = 0
    % The identity matrix is assigned to tensor diffussion matrix to avoid D=zeros   
    [szx szy szz szH] = size(T);
    sumvecT = zeros([szx*szy*szz, 1]);
    counter_sumvecT = 1;
    anisotropic_cofficient = 0.8;
    iosotropic_vec = [1; 0; 0; 1; 0; 1];
    for i = 1 : szx
        for j = 1 : szy
            for k = 1 : szz
                d11 = T(i,j,k,1);
                d12 = T(i,j,k,2);
                d13 = T(i,j,k,3);
                d22 = T(i,j,k,4);
                d23 = T(i,j,k,5);
                d33 = T(i,j,k,6);
                T_vec = squeeze(T(i,j,k,:));
                hessianmat = hessianvaluetomat(T_vec);
                temp_sum  = abs(d11) + abs(d12) + abs(d13) + abs(d22) + abs(d23) + abs(d33);
                det_hessianmat = det(hessianmat);
                if ((temp_sum == 0) || (det_hessianmat == 0))
                    T(i,j,k,1) = 1; T(i,j,k,4) = 1; T(i,j,k,6) = 1;
                else
                    T_vec = T_vec / det_hessianmat^(1/3);
                    T(i,j,k,:) = anisotropic_cofficient * iosotropic_vec + (1 - anisotropic_cofficient) *  T_vec;  
                end
%                 T(i,j,k,1) = 1; T(i,j,k,4) = 1; T(i,j,k,6) = 1;
%                 T(i,j,k,2) = 0; T(i,j,k,3) = 0; T(i,j,k,5) = 0;
            end
        end
    end
    boundary_value = zeros(szI(1),szI(2),szI(3));
    object = zeros(szI(1),szI(2),szI(3));
    object(SourcePoint(1),SourcePoint(2),SourcePoint(3)) = 1;
    volDim = [1,1,1];
    T_map = mxAnisoDistanceTransform(object, T, boundary_value, SpeedImage, volDim);
end