function [p, q, r] = afmMatrixInverse3by3(a, b, c)

	matA = [a(1), a(2), a(3); b(1), b(2), b(3); c(1), c(2), c(3)];
	detAmat = det(matA);
	if (detAmat == 0)
		fprintf('The MATRIX is singular, please check your Matrix');
		return;
	end
	% p(1) = (b(2)*c(3)-b(3)*c(2))/detAmat;
 %  	p(2) = (a(3)*c(2)-a(2)*c(3))/detAmat;
 %  	p(3) = (a(2)*b(3)-a(3)*b(2))/detAmat;
 %  	q(1) = (b(3)*c(1)-b(1)*c(3))/detAmat;
 %  	q(2) = (a(1)*c(3)-a(3)*c(1))/detAmat;
 %  	q(3) = (a(3)*b(1)-a(1)*b(3))/detAmat;
 %  	r(1) = (b(1)*c(2)-b(2)*c(1))/detAmat;
 %  	r(2) = (a(2)*c(1)-a(1)*c(2))/detAmat;
 %  	r(3) = (a(1)*b(2)-a(2)*b(1))/detAmat;
  	invmatA = inv(matA);
  	p = invmatA(1,:);
  	q = invmatA(2,:);
  	r = invmatA(3,:);
	% The following line is used to check the determinant calculation.  
	% fprintf('determinant calculated by matlab function is : %d\n', detAmat);
	% detA = a(1)*(b(2)*c(3) - b(3)*c(2)) - a(2)*(b(1)*c(3) - b(3)*c(1)) + a(3)*(b(1)*c(2) - b(2)*c(1));
	% fprintf('determinant calculated by afm method is : %d\n', detA);


end
% 	void MatrixInverse3by3( float* a, float* b, float* c, float* p, float* q, float* r )
% {
%   // FIRST THREE VECTORS ARE THE ROW VECTORS OF MATRIX A
%   // SECOND THREE VECTORS ARE THE ROW VECTORS OF MATRIX INV(A)
%   float detA = a[0]*(b[1]*c[2] - b[2]*c[1]) - a[1]*(b[0]*c[2] - b[2]*c[0]) + a[2]*(b[0]*c[1] - b[1]*c[0]);
%   if(detA == 0)
%     {
%       std::cout<<"["<<a[0]<<","<<a[1]<<","<<a[2]<<"]"<<std::endl;
%       std::cout<<"["<<b[0]<<","<<b[1]<<","<<b[2]<<"]"<<std::endl;
%       std::cout<<"["<<c[0]<<","<<c[1]<<","<<c[2]<<"]"<<std::endl;

%       std::cerr<<"Matrix singular!"<<std::endl;
%       exit(1);
%     }
%   p[0] = (b[1]*c[2]-b[2]*c[1])/detA;
%   p[1] = (a[2]*c[1]-a[1]*c[2])/detA;
%   p[2] = (a[1]*b[2]-a[2]*b[1])/detA;
%   q[0] = (b[2]*c[0]-b[0]*c[2])/detA;
%   q[1] = (a[0]*c[2]-a[2]*c[0])/detA;
%   q[2] = (a[2]*b[0]-a[0]*b[2])/detA;
%   r[0] = (b[0]*c[1]-b[1]*c[0])/detA;
%   r[1] = (a[1]*c[0]-a[0]*c[1])/detA;
%   r[2] = (a[0]*b[1]-a[1]*b[0])/detA;
% }