clear all; close all;
fprintf('Prepare for the original image : I binary image : bI\n');
% load('zebraI.mat');

% The first input variable I is inside op1resample.mat 
load('op1resample.mat');

% The second input variable is plot
plot_value = false;

% The third input variable is percentage 
percentage = 0.98;

% The fourth input vairable is rewire
rewire =  false;

% The fifth input vairable is gap 
gap = 10;

% The sixth input variable is ax
ax_value = false;

% The seventh input variable is dumpcheck
dumpcheck = true;

% The eighth input variable is connectrate
connectrate = 1.2;

% The ninth input variable is branchlen
branchlen = 8;

% The tenth input variable is branchlen
somagrowthcheck = false;

% The twelfth input variable is branchlen
cleanercheck = false;

threshold = 89.5;
bI = I > threshold;
 
fprintf('Prepare for the speed image.\n');

disp('Distance transform');
notbI = not(bI);
bdist = bwdist(notbI, 'Quasi-Euclidean');
bdist = bdist .* double(bI);
bdist = double(bdist);

disp('Looking for the source point...')
if somagrowthcheck
    SourcePoint = [soma.x; soma.y; soma.z];
    somaidx = find(soma.I == 1);
    [somax, somay, somaz] = ind2sub(size(soma.I), somaidx);
    % Find the soma radius
    d = pdist2([somax, somay, somaz], [soma.x, soma.y, soma.z]);
    maxD = max(d);
else
    [SourcePoint, maxD] = maxDistancePoint(bdist, I, true);
end

disp('Make the speed image...');
SpeedImage=(bdist/maxD).^4;
clear bdist;    
SpeedImage(SpeedImage==0) = 1e-10;

% 
% figure(1), imagesc(squeeze(max(SpeedImage,[],3))), title('speed');
% 
tic;
disp('marching...');
T = msfm(SpeedImage, SourcePoint, false, false);
% T = mxAnisoDistanceTransform(object, T, boundary, SpeedImage, volDim);
toc
disp('Saving the time crossing map for testing...');
save('T_script.mat','T')

% 
% Calculate gradient of DistanceMap
disp('Calculating gradient...')
grad = distgradient(T);

    % if plot
    %     axes(ax);
    % end
    S = {};
    B = zeros(size(T));
    if somagrowthcheck
        B = B | (soma.I>0.5);
    end
    lconfidence = [];
    if plot_value
        [x,y,z] = sphere;
        plot3(x + SourcePoint(2), y + SourcePoint(1), z + SourcePoint(3), 'ro');
    end

    unconnectedBranches = {};
    printcount = 0;
    printn = 0;

    while(true)

        StartPoint = maxDistancePoint(T, I, true);
        
        if plot_value
            plot3(x + StartPoint(2), y + StartPoint(1), z + StartPoint(3), 'ro');
        end

        if T(StartPoint(1), StartPoint(2), StartPoint(3)) == 0 || I(StartPoint(1), StartPoint(2), StartPoint(3)) == 0
            break;
        end

        [l, dump, merged, somamerged] = shortestpath2(T, grad, I, tree, StartPoint, SourcePoint, 1, 'rk4', gap);

        if size(l, 1) == 0
            l = StartPoint'; % Make sure the start point will be erased
        end
        % Get radius of each point from distance transform
        radius = zeros(size(l, 1), 1);
        parfor r = 1 : size(l, 1)
            radius(r) = getradius(I, l(r, 1), l(r, 2), l(r, 3));
        end
        radius(radius < 1) = 1;
        assert(size(l, 1) == size(radius, 1));

        [covermask, centremask] = binarysphere3d(size(T), l, radius);
        % Remove the traced path from the timemap
        if cleanercheck & size(l, 1) > branchlen
            covermask = augmask(covermask, I, l, radius);
        end

        % covermask(StartPoint(1), StartPoint(2), StartPoint(3)) = 3; % Why? Double check if it is nessensary - SQ

        T(covermask) = -1;
        T(centremask) = -3;

        % if cleanercheck
     %        T(wash==1) = -1;
     %    end

        % Add l to the tree
        if ~((dump) && dumpcheck) 
            [tree, newtree, conf, unconnected] = addbranch2tree(tree, l, merged, connectrate, radius, I, branchlen, plot, somamerged);
            lconfidence = [lconfidence, conf];
        end

        B = B | covermask;

        percent = sum(B(:) & I(:)) / sum(I(:));
%         if plot
%             axes(ax);
%         end
        printn = printn + 1;
        if printn > 1
            fprintf(1, repmat('\b',1,printcount));
            printcount = fprintf('Tracing percent: %f%%\n', percent*100);
        end
        if percent >= percentage
            disp('Coverage reached end tracing...')
            break;
        end


    end
% 
% % meanconf = mean(lconfidence);
% % 
% % if cleanercheck
% %     disp('Fixing topology')
% %     tree = fixtopology(tree);
% % end 
% % tree = prunetree(tree, branchlen);