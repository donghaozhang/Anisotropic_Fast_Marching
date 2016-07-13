function limitervalue = afmlimiter(x)
% float limiter( float x )
% {
%   float eps = 1e-7;
%   if(fabs(x)>eps)
%     return(x);
%   else
%     return(0);
% }
	afmeps = 1.0000e-07;
	if (afmabs_float(x) > afmeps)
		limitervalue = x;
	else
		limitervalue = 0;
	end
end

