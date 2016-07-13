function minvalue = afmmin(a, b) 
% original afm min code
% float min(float a, float b)
% {
%   if(a>=b)
%     return b;
%   else
%     return a;

% }
	if (a>=b)
		minvalue = b;
	else
		minvalue = a;
	end
end