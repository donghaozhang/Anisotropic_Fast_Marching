%
%Author  : Ender Konukoglu
%Email   : ender.konukoglu@gmail.com
%Article : @inproceedings{konukoglu2007recursive,
%          title={A recursive anisotropic fast marching approach to reaction diffusion equation: Application to tumor growth modeling},
%          author={Konukoglu, Ender and Sermesant, Maxime and Clatz, Olivier and Peyrat, Jean-Marc and Delingette, Herve and Ayache, Nicholas},
%          booktitle={Information processing in medical imaging},
%          pages={687--699},
%          year={2007},
%          organization={Springer}
%					}
%Date    : June 21, 2013.

%
 
function y=gaussian(x,m,sigma)
% m is the mean vector
% sigma is the covariance matrix
% x is the observation
% y is the resulting probability
y=1/((2*pi)^(length(x)/2)*sqrt(det(sigma)))*exp(-1/2*(x-m)'*inv(sigma)*(x-m));
