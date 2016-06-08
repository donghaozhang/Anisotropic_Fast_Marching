load('C:\Users\donghao\Desktop\Rivulet-Neuron-Tracing-Toolbox\SpeedImage_rivulet.mat')
rivulet_script = load('C:\Users\donghao\Desktop\AnisotropicFastMarchingexp - Copy\SpeedImage_script.mat');
rivulet_script = rivulet_script.SpeedImage;
SpeedImage_vec = SpeedImage(:); 
SpeedImage_script_vec = rivulet_script(:);
test_vec = SpeedImage_vec == SpeedImage_script_vec;
sum(test_vec(:))
figure(1),imagesc(squeeze(max(SpeedImage,[],3))), title('original rivulet');
figure(2),imagesc(squeeze(max(rivulet_script,[],3))), title('script');