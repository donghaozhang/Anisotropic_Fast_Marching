function subscriptpt = afmind2sub(afmsize, afmindex)
	[subscriptpt.x subscriptpt.y subscriptpt.z] = ind2sub(afmsize, afmindex);
    % [I,J] = ind2sub(SIZ,IND) returns the arrays I and J containing the
    % equivalent row and column subscripts corresponding to the index
    % matrix IND for a matrix of size SIZ. 
% The following code is afm ind2sub
%     void ind2sub(double* x, double* y, double* z, double ind, int* Size)
% {
%   float fracPart;
%   fracPart = modf( ind/(double)(Size[0]*Size[1]), z );
%   fracPart = fracPart*Size[1];
%   fracPart = modf( fracPart, y );
%   *x = fracPart*Size[0];
% }
end