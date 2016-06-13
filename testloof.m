clear opts; 
opts.useabsolute = 0; 
opts.responsetype = 2; 
opts.normalizationtype = 1;
load('zebraI.mat');
for i = 1:10
	radius = [1:i];
	hessian_out = oof3response(double(I), radius, opts);
	prefix_oof = 'mat\oof_r_';
	suffix_oof = '.mat';
	save([prefix_oof num2str(numel(radius)) suffix_oof],'hessian_out');
end