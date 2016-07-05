function [outputfeature_11, outputfeature_12, outputfeature_13, outputfeature_22, outputfeature_23, outputfeature_33] = oofhessian(image, radii)
	    tempoof = zeros(size(image, 1), size(image, 2), size(image, 3), numel(radii));
	    leig1 = zeros(size(image, 1), size(image, 2), size(image, 3), numel(radii));
	    leig2 = zeros(size(image, 1), size(image, 2), size(image, 3), numel(radii));
	    leig3 = zeros(size(image, 1), size(image, 2), size(image, 3), numel(radii));

	    % marginwidth is used to define the image border which is cropped to 
	    % circumvent the FFT wrap-around artifacts and save memory during 
	    % calculation.
	    % It is disabled by default. Enable it by uncommenting the second line
	    % below this description.
	    marginwidth = [0 0 0];
	    % marginwidth = [ceil((max(radii)+sigma*3)/pixelspacing(1)) ceil((max(radii)+sigma*3)/pixelspacing(2)) ceil((max(radii)+sigma*3)/pixelspacing(3)) ];

	    output=image(marginwidth(1)+1:end-marginwidth(1), marginwidth(2)+1:end-marginwidth(2), marginwidth(3)+1:end-marginwidth(3))*0; 
	    
	    % Default options
	    rtype = 0;
	    etype = 1;
	    ntype = 1;
	    pixelspacing=[1 1 1];
	    sigma= min(pixelspacing);
	   
	    if exist('options', 'var')~=0
	        if isfield(options, 'spacing')~=0
	            pixelspacing = options.spacing;
	            sigma= min(pixelspacing);
	        end
	        if isfield(options, 'responsetype')~=0
	            rtype = options.responsetype;
	        end
	        if isfield(options, 'normalizationtype')~=0
	            ntype = options.normalizationtype;
	        end
	        if isfield(options, 'sigma')~=0
	            sigma = options.sigma;
	        end
	        if isfield(options, 'useabsolute')~=0
	            etype = options.useabosolute;
	        end
	        if ((min(radii)<sigma) & ntype>0)
	            disp('Sigma must be >= minimum range to enable the advanced normalization. The current setting falls back to options.normalizationtype=0, because of the undersize sigma.');
	            ntype = 0;
	        end
	    end
	    
	    imgfft = fftn(image);
	    
	    %Obtaining the Fourier coordinate
	    [x,y,z] = ifftshiftedcoormatrix([size(image,1) size(image,2) size(image,3)]);
	    
	    
	    % Casting the Fourier coordiantes to be of the same type as image
	    x=x+image(1)*0;
	    y=y+image(1)*0;
	    z=z+image(1)*0;
	    % End of the type casting
	    x=x/size(image,1)/pixelspacing(1);
	    y=y/size(image,2)/pixelspacing(2);
	    z=z/size(image,3)/pixelspacing(3);        
	    radius=realsqrt(x.^2+y.^2+z.^2)+1e-12;
	    
	    % Save memory by clearing x y z. Although obtained from different
	    % functions, x y z are equivalent to:
	    % x = ifftshiftedcoordinate(size(image), 1, pixelspacing)
	    % y = ifftshiftedcoordinate(size(image), 2, pixelspacing)
	    % z = ifftshiftedcoordinate(size(image), 3, pixelspacing)
	    % If main memory (or GPU memory) has enough memory to buffer the 
	    % entire x,y,z, comment the following clear command and replace the 
	    % equivalent bufferred variables inside the following for-loop. It 
	    % gives around 20% speed up.
	    clear x y z 

	    for i=1:length(radii)
	        fprintf('Working on radii %f\n', radii(i));
	        normalization = 4/3*pi*radii(i)^3/(besselj(1.5, 2*pi*radii(i)*1e-12)/(1e-12)^(3/2)) /radii(i)^2 * ((radii(i)/sqrt(2*radii(i)*sigma-sigma*sigma))^ntype);
	        
	        besseljBuffer = normalization * exp((-(sigma)^2)*2*pi*pi* (radius.^2))./(radius.^(3/2));
	        besseljBuffer = ( sin(2*pi*radii(i)*radius)./(2*pi*radii(i)*radius) - cos(2*pi*radii(i)*radius)) .* besseljBuffer.*sqrt(1/pi/pi/radii(i)./radius) ;

	        % clear radius
	        besseljBuffer=besseljBuffer.*imgfft;

	% There are 6 3D IFFT performed at each radius. Here we use in-place FFT to
	% save memory, although the code looks clumpsy.
	% If you are using Cuda-enabled GPU acceleration or you have large enough 
	% memory to use out-of-place FFT, uncomment the following 6 lines, and
	% comment the inplace FFT codes. It gives about 20%-40% overall speed up.
	%          outputfeature_11 = freqOp(real(ifftn(x.*x.*besseljBuffer)), marginwidth);
	%          outputfeature_12 = freqOp(real(ifftn(x.*y.*besseljBuffer)), marginwidth); 
	%          outputfeature_13 = freqOp(real(ifftn(x.*z.*besseljBuffer)), marginwidth);
	% % 
	%          outputfeature_22 = freqOp(real(ifftn(y.*y.*besseljBuffer)), marginwidth); 
	%          outputfeature_23 = freqOp(real(ifftn(y.*z.*besseljBuffer)), marginwidth);
	% % 
	%          outputfeature_33 = freqOp(real(ifftn(z.*z.*besseljBuffer)), marginwidth); 

	% Inplace FFT
	        buffer=ifftshiftedcoordinate(size(image), 1, pixelspacing).^2.* ............x.*x.*  .....    
	                                     besseljBuffer;
	        buffer=ifft(buffer, [], 1);buffer=ifft(buffer, [], 2);buffer=ifft(buffer, [], 3, 'symmetric');
	        buffer = freqOp(buffer, marginwidth); outputfeature_11 = buffer;
	        clear buffer;
	        buffer=ifftshiftedcoordinate(size(image), 1, pixelspacing).*ifftshiftedcoordinate(size(image), 2, pixelspacing).* ........x.*y.*  .....    
	                                     besseljBuffer;
	        buffer=ifft(buffer, [], 1);buffer=ifft(buffer, [], 2);buffer=ifft(buffer, [], 3, 'symmetric');                                 
	        buffer = freqOp(buffer, marginwidth); outputfeature_12 = buffer;
	        clear buffer;        
	        buffer=ifftshiftedcoordinate(size(image), 1, pixelspacing).*ifftshiftedcoordinate(size(image), 3, pixelspacing).* ........x.*z.*  .....    
	                                     besseljBuffer;                                 
	        buffer=ifft(buffer, [], 1);buffer=ifft(buffer, [], 2);buffer=ifft(buffer, [], 3, 'symmetric');                                 
	        buffer = freqOp(buffer, marginwidth); outputfeature_13 = buffer;
	        clear buffer;        

	        buffer=ifftshiftedcoordinate(size(image), 2, pixelspacing).^2.* .........*y.*y  .....    
	                                     besseljBuffer;
	        buffer=ifft(buffer, [], 1);buffer=ifft(buffer, [], 2);buffer=ifft(buffer, [], 3, 'symmetric');                                 
	        buffer = freqOp(buffer, marginwidth); outputfeature_22 = buffer;
	        clear buffer;        
	        buffer=ifftshiftedcoordinate(size(image), 2, pixelspacing).*ifftshiftedcoordinate(size(image), 3, pixelspacing).* ........y.*z.*  .....    
	                                     besseljBuffer;
	        buffer=ifft(buffer, [], 1);buffer=ifft(buffer, [], 2);buffer=ifft(buffer, [], 3, 'symmetric');                                 
	        buffer = freqOp(buffer, marginwidth); outputfeature_23 = buffer;
	        clear buffer;            

	        buffer=ifftshiftedcoordinate(size(image), 3, pixelspacing).^2.* ........ z.*z.*  .....    
	                                     besseljBuffer;
	        buffer=ifft(buffer, [], 1);buffer=ifft(buffer, [], 2);buffer=ifft(buffer, [], 3, 'symmetric');        
	        buffer = freqOp(buffer, marginwidth); outputfeature_33 = buffer;      
	    end
end