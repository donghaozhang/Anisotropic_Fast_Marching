function [R, deltalogic] = afmfind_roots_indic(a, b, c)
% bool find_roots_indic(float* R,float a,float b, float c)
% {


%   if(limiter(b*b - 4*a*c) < 0)
%     {
%       R[0] = 0;
%       R[1] = 0;
%       return(true);
%     }
%   else
%     {
      
%       R[0] = -b/(2*a) + sqrt( limiter(b*b - 4*a*c) )/(2*a);
%       R[1] = -b/(2*a) - sqrt( limiter(b*b - 4*a*c) )/(2*a);
%       return(false);
%     }


% }
% matlab way
% p = [3 -2 -4];
% r = roots(p)
  if (afmlimiter(b*b - 4*a*c) < 0)
    R(1) = 0;
    R(2) = 0;
    deltalogic = true;
  else
    R(1) = -b / (2 * a) + sqrt(afmlimiter(b * b - 4 * a * c))/(2*a);
    R(2) = -b / (2 * a) - sqrt(afmlimiter(b * b - 4 * a * c))/(2*a);
    deltalogic = false; 
  end
end
