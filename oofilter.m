% clear all; 
close all;
% fprintf('Prepare for the original image : I binary image : bI\n');
clear opts; 
opts.useabsolute = 0; 
opts.responsetype = 3; 
opts.normalizationtype = 1;
% The following line is coded for zebrafish 
% load('midzebraI.mat');
% load mouseRGC image mat 
% load('miccaidata\experimentmouseRGC\mouseRGCresampled.mat');
% I = permute(I, [2 1 3]); 
% load('op1resample.mat');
radius = 1 : 3; 
hessian_out = oof3response(double(I), radius, opts);
% midzebra oof response image threshold
% thresholdoof = 0.42; 
% the threshold of mouseRGC
% thresholdoof = 0.18;
% the threhold of janelia fly
thresholdoof = 2.4;
B = hessian_out > thresholdoof; 
% clear opts;
clear radius;
clear x;
clear y;
clear z;
clear thresholdoof;
thresI = I > threshold;
% I try to combine the binary image generated the simple threshold and optimnal oriented flux filter
combineB = B | thresI;
% combineB = B;
[x y z] = ind2sub(size(combineB), find(combineB));
plot3(y, x, z, '.', 'Color', [0.7 0.7 1]);  view(3); axis equal;
I = combineB;
clear B;
clear radius;
clear hessian_out;
clear thresI;
clear opts;
% [Iout,whatScale,Voutx,Vouty,Voutz]=FrangiFilter3D(I_original);