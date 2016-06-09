clear opts; 
opts.useabsolute = 0; 
opts.responsetype = 1; 
opts.normalizationtype = 0;
load('zebraI.mat');
for i = 1 : 10
	 radius = [1:i];
	 oof_result = learnoof(double(I), radius, opts);
	 prefix_oof = 'oof_r_';
	 suffix_oof = '.mat';
	save([prefix_oof num2str(numel(radius)) suffix_oof],'oof_result');
end