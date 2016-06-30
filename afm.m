function T_map = afm(I, threshold, foreground_speed_coefficient, speedastensorflag, oofhmflag, anisotropic_cofficient)
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
%     figure(1),imagesc(squeeze(max(SpeedImage,[],3))'), title('speed iamge xy projection rivulet');
%     figure(2),imagesc(squeeze(max(SpeedImage,[],2))), title('speed iamge xy projection original rivulet');
    % Original 
    % SpeedImage=(bdist/maxD).^4;
    % clear bdist;
    % SpeedImage(SpeedImage==0) = 1e-10;
    szI = size(I);
    if (~oofhmflag)
        sigma_value = 0.7;
        if speedastensorflag
            [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(double(bdist), sigma_value);
        else
            [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(double(I), sigma_value);
            fprintf('The hessian matrix is derived from frangi.\n');
        end

        T = zeros(szI(1),szI(2),szI(3),6);

        eps = 1e-5;
        scale_para = 3;
        T(:,:,:,1) = Dxx * scale_para;
        T(:,:,:,2) = Dxy * scale_para;
        T(:,:,:,3) = Dxz * scale_para;
        T(:,:,:,4) = Dyy * scale_para;
        T(:,:,:,5) = Dyz * scale_para;
        T(:,:,:,6) = Dzz * scale_para;
        save('mat\hmatvess.mat','T');
    elseif oofhmflag
        clear opts; 
        opts.useabsolute = 0; 
        opts.responsetype = 1; 
        opts.normalizationtype = 0;
        radius = [1:7];
        T = oof_hessian(double(I), radius, opts);
        fprintf('The hessian matrix is derived from optimal oriented flux.\n');
        save('mat\hmatoof.mat','T');
    end
            
    % why I do the following code is make sure that Dxx = 1; Dyy = 1; Dzz = 1;
    % Dxy == Dyx = 0 Dxz == Dzx = 0
    % The identity matrix is assigned to tensor diffussion matrix to avoid D=zeros   
    [szx szy szz szH] = size(T);
    sumvecT = zeros([szx*szy*szz, 1]);
    counter_sumvecT = 1;
    % anisotropic_cofficient = 0.95;
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
                    % T_vec = T_vec / ((abs(det_hessianmat))^(1/3)) / 20;
                    T_vec = T_vec / norm(T_vec) * 10;
                    % T(i,j,k,:) = anisotropic_cofficient * iosotropic_vec + (1 - anisotropic_cofficient) *  T_vec;  
                    T(i,j,k,:) = (1 - anisotropic_cofficient) *  T_vec;  
                end
%                 T(i,j,k,1) = 1; T(i,j,k,4) = 1; T(i,j,k,6) = 1;
%                 T(i,j,k,2) = 0; T(i,j,k,3) = 0; T(i,j,k,5) = 0;
            end
        end
    end
    save('mat\diffussion.mat','T');
    boundary_value = zeros(szI(1),szI(2),szI(3));
    object = zeros(szI(1),szI(2),szI(3));
    object(SourcePoint(1),SourcePoint(2),SourcePoint(3)) = 1;
    volDim = [1,1,1];
    T_map = mxAnisoDistanceTransform(object, T, boundary_value, SpeedImage, volDim);
end