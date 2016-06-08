load('C:\Users\donghao\Desktop\Rivulet-Neuron-Tracing-Toolbox\bdist_rivulet.mat')
rivulet_script = load('C:\Users\donghao\Desktop\AnisotropicFastMarchingexp - Copy\bdist_script.mat');
rivulet_script = rivulet_script.bdist;
bdist_vec = bdist(:); 
bdist_script_vec = rivulet_script(:);
test_vec = bdist_vec == bdist_script_vec;
sum(test_vec(:))
figure(1),imagesc(squeeze(max(bdist,[],3))), title('original rivulet');
figure(2),imagesc(squeeze(max(rivulet_script,[],3))), title('script');