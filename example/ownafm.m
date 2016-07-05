tic
clear all; close all;
fprintf('Prepare for the original image : I binary image : bI\n');
% load('zebraI.mat');
load('op1resample.mat');
threshold = 20;
bI = I > threshold;

fprintf('Prepare for the speed image.\n');

disp('Distance transform');
notbI = not(I>threshold);
bdist = bwdist(notbI, 'Quasi-Euclidean');
bdist = bdist .* double(I);
bdist = double(bdist);
[SourcePoint, maxD] = maxDistancePoint(bdist, I, true);
% Speical treatment for anisotropic fast marching
SpeedImage= bdist/maxD * 10;
SpeedImage(SpeedImage==0) = 0.1;

% Original 
% SpeedImage=(bdist/maxD).^4;
% clear bdist;
% SpeedImage(SpeedImage==0) = 1e-10;

sigma_value = 0.42;
% SpeedImage = double(bI) * 2;

figure(1), imagesc(squeeze(max(SpeedImage,[],3))), title('speed');

[DX,DY,DZ] = gradient(double(I));
% imagesc(max(DX,[],3))
[Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(double(I), sigma_value);
szI = size(I);
T = zeros(szI(1),szI(2),szI(3),6);

eps = 1e-5;
T(:,:,:,1) = Dxx;
T(:,:,:,2) = Dxy;
T(:,:,:,3) = Dxz;
T(:,:,:,4) = Dyy;
T(:,:,:,5) = Dyz;
T(:,:,:,6) = Dzz;
% why I do the following code is make sure that Dxx = 1; Dyy = 1; Dzz = 1;
% Dxy == Dyx = 0 Dxz == Dzx = 0
% The identity matrix is assigned to tensor diffussion matrix to avoid D=zeros   
[szx szy szz szH] = size(T);
sumvecT = zeros([szx*szy*szz, 1]);
counter_sumvecT = 1;
for i = 1 : szx
    for j = 1 : szy
        for k = 1 : szz
            d11 = T(i,j,k,1);
            d12 = T(i,j,k,2);
            d13 = T(i,j,k,3);
            d22 = T(i,j,k,4);
            d23 = T(i,j,k,5);
            d33 = T(i,j,k,6);
            temp_sum  = d11 + d12 + d13 + d22 + d23 + d33;
            if (temp_sum == 0)
                T(i,j,k,1) = 1; T(i,j,k,4) = 1; T(i,j,k,6) = 1;
            end    
        end
    end
end
boundary = zeros(szI(1),szI(2),szI(3));
object = zeros(szI(1),szI(2),szI(3));
object(SourcePoint(1),SourcePoint(2),SourcePoint(3)) = 1;
volDim = [1,1,1];
% we use our own speed image instead of one provided by the author of
% anisotropic fast marching

tic;
% Distance = msfm(SpeedImage, SourcePoint, false, false);
Distance = mxAnisoDistanceTransform(object, T, boundary, SpeedImage, volDim);
toc

figure(2),imagesc(squeeze(max(Distance,[],3))), title('XY max projection');
d = zeros(size(Distance));
d(Distance > 0) = 1;
% d = d ./ double(bI);
figure(3),imagesc(squeeze(max(d,[], 3))), title('T tseg');

% % Calculate gradient of DistanceMap
% disp('Calculating gradient...')
% grad = distgradient(Distance);
% unconnectedBranches = {};
% printcount = 0;
% printn = 0;
% % The following line initialises the swc tree
% tree = [];
% dumpcheck = false;
% plot_value = false;
% gap = 10;
% ax = false;
% cleanercheck = false;
% percentage = 0.9;
% connectrate = 1.5;
% branchlen = 4;
% dump = false;
% lconfidence = [];
% S = {};
% B = zeros(size(Distance));
% % while(true)
% for i = 1 : 40000
%     StartPoint = maxDistancePoint(Distance, I, true);
%     fprintf('Stop Position one: \n');
%     if plot_value
% 	    plot3(x + StartPoint(2), y + StartPoint(1), z + StartPoint(3), 'ro');
%     end
%     fprintf('Stop Position two: \n');
%     if Distance(StartPoint(1), StartPoint(2), StartPoint(3)) == 0 || I(StartPoint(1), StartPoint(2), StartPoint(3)) == 0
%     	break;
%     end

%     [l, dump, merged, somamerged] = shortestpath2(Distance, grad, I, tree, StartPoint, SourcePoint, 1, 'rk4', gap);
%     if size(l, 1) == 0
%         l = StartPoint'; % Make sure the start point will be erased
%     end
%     % Get radius of each point from distance transform
%     radius = zeros(size(l, 1), 1);
%     parfor r = 1 : size(l, 1)
% 	    radius(r) = getradius(I, l(r, 1), l(r, 2), l(r, 3));
%     end
%     fprintf('Stop Position three: \n');
%     radius(radius < 1) = 1;
% 	assert(size(l, 1) == size(radius, 1));

%     [covermask, centremask] = binarysphere3d(size(Distance), l, radius);
%     % Remove the traced path from the timemap
%     if cleanercheck & size(l, 1) > branchlen
%         covermask = augmask(covermask, I, l, radius);
%     end

%     % covermask(StartPoint(1), StartPoint(2), StartPoint(3)) = 3; % Why? Double check if it is nessensary - SQ

%     Distance(covermask) = -1;
%     Distance(centremask) = -3;

%     % if cleanercheck
%  %        T(wash==1) = -1;
%  %    end

%     % Add l to the tree
%     if ~((dump) && dumpcheck) 
% 	    [tree, newtree, conf, unconnected] = addbranch2tree(tree, l, merged, connectrate, radius, I, branchlen, plot_value, somamerged);
%         lconfidence = [lconfidence, conf];
% 	end

%     B = B | covermask;

%     percent = sum(B(:) & I(:)) / sum(I(:));
%     if plot_value
%         axes(ax);
%     end
%     printn = printn + 1;
%     if printn > 1
%         fprintf(1, repmat('\b',1,printcount));
%         printcount = fprintf('Tracing percent: %f%%\n', percent*100);
%     end
%     if percent >= percentage
%     	disp('Coverage reached end tracing...')
%     	break;
%     end

% end

% meanconf = mean(lconfidence);
% 
% if cleanercheck
%     disp('Fixing topology')
%     tree = fixtopology(tree);
% end 
% tree = prunetree(tree, branchlen);

% %% Form now on, I am going to use oof to extract diffusion tensor matrix
% clear all;close all;
% fprintf('Prepare for the original image : I binary image : bI\n');
% % load('zebraI.mat');
% load('op1resample.mat');
% image = double(I);
% radii = [1:0.5:4];
%     tempoof = zeros(size(image, 1), size(image, 2), size(image, 3), numel(radii));
%     leig1 = zeros(size(image, 1), size(image, 2), size(image, 3), numel(radii));
%     leig2 = zeros(size(image, 1), size(image, 2), size(image, 3), numel(radii));
%     leig3 = zeros(size(image, 1), size(image, 2), size(image, 3), numel(radii));
% 
%     % marginwidth is used to define the image border which is cropped to 
%     % circumvent the FFT wrap-around artifacts and save memory during 
%     % calculation.
%     % It is disabled by default. Enable it by uncommenting the second line
%     % below this description.
%     marginwidth = [0 0 0];
%     % marginwidth = [ceil((max(radii)+sigma*3)/pixelspacing(1)) ceil((max(radii)+sigma*3)/pixelspacing(2)) ceil((max(radii)+sigma*3)/pixelspacing(3)) ];
% 
%     output=image(marginwidth(1)+1:end-marginwidth(1), marginwidth(2)+1:end-marginwidth(2), marginwidth(3)+1:end-marginwidth(3))*0; 
%     
%     % Default options
%     rtype = 0;
%     etype = 1;
%     ntype = 1;
%     pixelspacing=[1 1 1];
%     sigma= min(pixelspacing);
%    
%     if exist('options', 'var')~=0
%         if isfield(options, 'spacing')~=0
%             pixelspacing = options.spacing;
%             sigma= min(pixelspacing);
%         end
%         if isfield(options, 'responsetype')~=0
%             rtype = options.responsetype;
%         end
%         if isfield(options, 'normalizationtype')~=0
%             ntype = options.normalizationtype;
%         end
%         if isfield(options, 'sigma')~=0
%             sigma = options.sigma;
%         end
%         if isfield(options, 'useabsolute')~=0
%             etype = options.useabosolute;
%         end
%         if ((min(radii)<sigma) & ntype>0)
%             disp('Sigma must be >= minimum range to enable the advanced normalization. The current setting falls back to options.normalizationtype=0, because of the undersize sigma.');
%             ntype = 0;
%         end
%     end
%     
%     imgfft = fftn(image);
%     
%     %Obtaining the Fourier coordinate
%     [x,y,z] = ifftshiftedcoormatrix([size(image,1) size(image,2) size(image,3)]);
%     
%     
%     % Casting the Fourier coordiantes to be of the same type as image
%     x=x+image(1)*0;
%     y=y+image(1)*0;
%     z=z+image(1)*0;
%     % End of the type casting
%     x=x/size(image,1)/pixelspacing(1);
%     y=y/size(image,2)/pixelspacing(2);
%     z=z/size(image,3)/pixelspacing(3);        
%     radius=realsqrt(x.^2+y.^2+z.^2)+1e-12;
%     
%     % Save memory by clearing x y z. Although obtained from different
%     % functions, x y z are equivalent to:
%     % x = ifftshiftedcoordinate(size(image), 1, pixelspacing)
%     % y = ifftshiftedcoordinate(size(image), 2, pixelspacing)
%     % z = ifftshiftedcoordinate(size(image), 3, pixelspacing)
%     % If main memory (or GPU memory) has enough memory to buffer the 
%     % entire x,y,z, comment the following clear command and replace the 
%     % equivalent bufferred variables inside the following for-loop. It 
%     % gives around 20% speed up.
%     clear x y z 
% 
%     for i=1:length(radii)
%         fprintf('Working on radii %f\n', radii(i));
%         normalization = 4/3*pi*radii(i)^3/(besselj(1.5, 2*pi*radii(i)*1e-12)/(1e-12)^(3/2)) /radii(i)^2 * ((radii(i)/sqrt(2*radii(i)*sigma-sigma*sigma))^ntype);
%         
%         besseljBuffer = normalization * exp((-(sigma)^2)*2*pi*pi* (radius.^2))./(radius.^(3/2));
%         besseljBuffer = ( sin(2*pi*radii(i)*radius)./(2*pi*radii(i)*radius) - cos(2*pi*radii(i)*radius)) .* besseljBuffer.*sqrt(1/pi/pi/radii(i)./radius) ;
% 
%         % clear radius
%         besseljBuffer=besseljBuffer.*imgfft;
% 
% % There are 6 3D IFFT performed at each radius. Here we use in-place FFT to
% % save memory, although the code looks clumpsy.
% % If you are using Cuda-enabled GPU acceleration or you have large enough 
% % memory to use out-of-place FFT, uncomment the following 6 lines, and
% % comment the inplace FFT codes. It gives about 20%-40% overall speed up.
% %          outputfeature_11 = freqOp(real(ifftn(x.*x.*besseljBuffer)), marginwidth);
% %          outputfeature_12 = freqOp(real(ifftn(x.*y.*besseljBuffer)), marginwidth); 
% %          outputfeature_13 = freqOp(real(ifftn(x.*z.*besseljBuffer)), marginwidth);
% % % 
% %          outputfeature_22 = freqOp(real(ifftn(y.*y.*besseljBuffer)), marginwidth); 
% %          outputfeature_23 = freqOp(real(ifftn(y.*z.*besseljBuffer)), marginwidth);
% % % 
% %          outputfeature_33 = freqOp(real(ifftn(z.*z.*besseljBuffer)), marginwidth); 
% 
% % Inplace FFT
%         buffer=ifftshiftedcoordinate(size(image), 1, pixelspacing).^2.* ............x.*x.*  .....    
%                                      besseljBuffer;
%         buffer=ifft(buffer, [], 1);buffer=ifft(buffer, [], 2);buffer=ifft(buffer, [], 3, 'symmetric');
%         buffer = freqOp(buffer, marginwidth); outputfeature_11 = buffer;
%         clear buffer;
%         buffer=ifftshiftedcoordinate(size(image), 1, pixelspacing).*ifftshiftedcoordinate(size(image), 2, pixelspacing).* ........x.*y.*  .....    
%                                      besseljBuffer;
%         buffer=ifft(buffer, [], 1);buffer=ifft(buffer, [], 2);buffer=ifft(buffer, [], 3, 'symmetric');                                 
%         buffer = freqOp(buffer, marginwidth); outputfeature_12 = buffer;
%         clear buffer;        
%         buffer=ifftshiftedcoordinate(size(image), 1, pixelspacing).*ifftshiftedcoordinate(size(image), 3, pixelspacing).* ........x.*z.*  .....    
%                                      besseljBuffer;                                 
%         buffer=ifft(buffer, [], 1);buffer=ifft(buffer, [], 2);buffer=ifft(buffer, [], 3, 'symmetric');                                 
%         buffer = freqOp(buffer, marginwidth); outputfeature_13 = buffer;
%         clear buffer;        
% 
%         buffer=ifftshiftedcoordinate(size(image), 2, pixelspacing).^2.* .........*y.*y  .....    
%                                      besseljBuffer;
%         buffer=ifft(buffer, [], 1);buffer=ifft(buffer, [], 2);buffer=ifft(buffer, [], 3, 'symmetric');                                 
%         buffer = freqOp(buffer, marginwidth); outputfeature_22 = buffer;
%         clear buffer;        
%         buffer=ifftshiftedcoordinate(size(image), 2, pixelspacing).*ifftshiftedcoordinate(size(image), 3, pixelspacing).* ........y.*z.*  .....    
%                                      besseljBuffer;
%         buffer=ifft(buffer, [], 1);buffer=ifft(buffer, [], 2);buffer=ifft(buffer, [], 3, 'symmetric');                                 
%         buffer = freqOp(buffer, marginwidth); outputfeature_23 = buffer;
%         clear buffer;            
% 
%         buffer=ifftshiftedcoordinate(size(image), 3, pixelspacing).^2.* ........ z.*z.*  .....    
%                                      besseljBuffer;
%         buffer=ifft(buffer, [], 1);buffer=ifft(buffer, [], 2);buffer=ifft(buffer, [], 3, 'symmetric');        
%         buffer = freqOp(buffer, marginwidth); outputfeature_33 = buffer;      
%     end
% threshold = 20;
% bI = I > threshold;
% 
% fprintf('Prepare for the speed image.\n');
% 
% disp('Distance transform');
% notbI = not(I>threshold);
% bdist = bwdist(notbI, 'Quasi-Euclidean');
% bdist = bdist .* double(I);
% bdist = double(bdist);
% [SourcePoint, maxD] = maxDistancePoint(bdist, I, true);
% SpeedImage= bdist/maxD * 10;
% SpeedImage(SpeedImage==0) = 0.1;
% 
% sigma_value = 0.42;
% % SpeedImage = double(bI) * 2;
% 
% figure(1), imagesc(squeeze(max(SpeedImage,[],3))), title('speed');
% 
% szI = size(I);
% T = zeros(szI(1),szI(2),szI(3),6);
% 
% eps = 1e-5;
% T(:,:,:,1) = outputfeature_11;
% T(:,:,:,2) = outputfeature_12;
% T(:,:,:,3) = outputfeature_13;
% T(:,:,:,4) = outputfeature_22;
% T(:,:,:,5) = outputfeature_23;
% T(:,:,:,6) = outputfeature_33;
% 
% boundary = zeros(szI(1),szI(2),szI(3));
% object = zeros(szI(1),szI(2),szI(3));
% object(SourcePoint(1),SourcePoint(2),SourcePoint(3)) = 1;
% volDim = [1,1,1];
% % we use our own speed image instead of one provided by the author of
% % anisotropic fast marching
% 
% tic;
% Distance = mxAnisoDistanceTransform(object, T, boundary, SpeedImage, volDim);
% toc
% 
% figure(2),imagesc(squeeze(max(Distance,[],3))), title('XY max projection');
% d = zeros(size(Distance));
% d(Distance > 0) = 1;
% % d = d ./ double(bI);
% figure(3),imagesc(squeeze(max(d,[], 3))), title('T tseg');