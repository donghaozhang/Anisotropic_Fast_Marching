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

function se = gaussian_structural_element_2D(covariance_matrix,S)
% creates 3D gaussian structural element with the given covariance matrix (3x3x3) and given size,


m = [ceil(S(1)/2) ceil(S(2)/2)]';


for j=1:S(1)
    for k=1:S(2)

            
            x = [j,k]';
            se(j,k)=gaussian(x,m,covariance_matrix);
            
        end;
    end;



