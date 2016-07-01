load('C:\Users\donghao\Desktop\Rivulet-Neuron-Tracing-Toolbox\T_rivulet.mat')
rivulet_script = load('C:\Users\donghao\Desktop\AnisotropicFastMarchingexp - Copy\T_script.mat');
rivulet_script = rivulet_script.T;
T_vec = T(:); 
T_script_vec = rivulet_script(:);
test_vec = T_vec == T_script_vec;
sum(test_vec(:))
figure(1),imagesc(squeeze(max(T,[],3))), title('original rivulet');
figure(2),imagesc(squeeze(max(rivulet_script,[],3))), title('script');