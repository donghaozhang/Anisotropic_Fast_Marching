%% radius 1 - 4
tmp_mat = load('C:\Users\donghao\Desktop\Anisotropic-Fast-Marching\mat\oof_r_1.mat');
oof_result = tmp_mat.hessian_out;
tmp_mat = squeeze(max(oof_result,[],3));
tmp_mat = permute(tmp_mat, [2,1]);
figure(1)
subplot(2,2,1)
imagesc(tmp_mat)
title('radius 1');
tmp_mat = load('C:\Users\donghao\Desktop\Anisotropic-Fast-Marching\mat\oof_r_2.mat');
oof_result = tmp_mat.hessian_out;
tmp_mat = squeeze(max(oof_result,[],3));
tmp_mat = permute(tmp_mat, [2,1]);
subplot(2,2,2)
imagesc(tmp_mat) 
title('radius 2');
tmp_mat = load('C:\Users\donghao\Desktop\Anisotropic-Fast-Marching\mat\oof_r_3.mat');
oof_result = tmp_mat.hessian_out;
tmp_mat = squeeze(max(oof_result,[],3));
tmp_mat = permute(tmp_mat, [2,1]);
subplot(2,2,3)
imagesc(tmp_mat) 
title('radius 3');
tmp_mat = load('C:\Users\donghao\Desktop\Anisotropic-Fast-Marching\mat\oof_r_4.mat');
oof_result = tmp_mat.hessian_out;
tmp_mat = squeeze(max(oof_result,[],3));
tmp_mat = permute(tmp_mat, [2,1]);
subplot(2,2,4)
imagesc(tmp_mat) 
title('radius 4');
%% radius from 5 -8
tmp_mat = load('C:\Users\donghao\Desktop\Anisotropic-Fast-Marching\mat\oof_r_5.mat');
oof_result = tmp_mat.hessian_out;
tmp_mat = squeeze(max(oof_result,[],3));
tmp_mat = permute(tmp_mat, [2,1]);
figure(2)
subplot(2,2,1)
imagesc(tmp_mat)
title('radius 5');
tmp_mat = load('C:\Users\donghao\Desktop\Anisotropic-Fast-Marching\mat\oof_r_6.mat');
oof_result = tmp_mat.hessian_out;
tmp_mat = squeeze(max(oof_result,[],3));
tmp_mat = permute(tmp_mat, [2,1]);
subplot(2,2,2)
imagesc(tmp_mat) 
title('radius 6');
tmp_mat = load('C:\Users\donghao\Desktop\Anisotropic-Fast-Marching\mat\oof_r_7.mat');
oof_result = tmp_mat.hessian_out;
tmp_mat = squeeze(max(oof_result,[],3));
tmp_mat = permute(tmp_mat, [2,1]);
subplot(2,2,3)
imagesc(tmp_mat) 
title('radius 7');
tmp_mat = load('C:\Users\donghao\Desktop\Anisotropic-Fast-Marching\mat\oof_r_8.mat');
oof_result = tmp_mat.hessian_out;
tmp_mat = squeeze(max(oof_result,[],3));
tmp_mat = permute(tmp_mat, [2,1]);
subplot(2,2,4)
imagesc(tmp_mat) 
title('radius 8');