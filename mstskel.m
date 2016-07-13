clear all;
close all;
threshold = 56;
% load('../../Anisotropic-Fast-Marching/zebraI.mat');
load('zebraI.mat');
testvol = I > threshold;
%load testvol
tic
skel = Skeleton3D(testvol);
[xpoint ypoint zpoint] = ind2sub(size(skel), find(skel));
% [mp(50, 1); z], [mp(50, 2); mp(:, 2)], [mp(50, 3); mp(:, 3)];
% x = rand(100,1)*100; y = rand(100,1)*100; z = zeros(100,1);
% tree = MST_tree(1, x, y,  z, .5, 50, [], [], 'none');
for edgei = 1 : size(xpoint, 1)
	for edgej= 1 : size(xpoint, 1)
            edgemap(edgei, edgej) = (xpoint(edgei) - xpoint(edgej)) * (xpoint(edgei) - xpoint(edgej))...
            	+ (ypoint(edgei) - ypoint(edgej)) * (ypoint(edgei) - ypoint(edgej)) ...
            	+ (zpoint(edgei) - zpoint(edgej)) * (zpoint(edgei) - zpoint(edgej));
    end
end
pindexlist = []
for i = 1 : size(xpoint, 1)
	[sortarray, sortedindex] = sort(edgemap(i,:));
	pindex = find( edgemap(i,:) == sortarray(2));
	curpindex = pindex(1);
	pindexlist(i) = curpindex;
	% swctree(i, 1) = i;
	% swctree(i, 2) = 2;
	% swctree(i, 3) = xpoint;
	% swctree(i, 4) = ypoint;
	% swctree(i, 5) = zpoint;
	% swctree(i, 6) = 1;
	% swctree(i, 7) =  curpindex;
end
swctree = zeros(size(xpoint, 1), 7);
for i = 1 : size(xpoint, 1)
	swctree(i, 1) = i;
	swctree(i, 2) = 2;
	swctree(i, 3) = xpoint(i);
	swctree(i, 4) = ypoint(i);
	swctree(i, 5) = zpoint(i);
	swctree(i, 6) = 1;
	swctree(i, 7) = pindexlist(i);
end
showswc(swctree);
outfilename = 'swc\mstskel.swc';
saveswc(swctree, outfilename);

