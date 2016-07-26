function groupvvaluemat = afmgroup_velocity(x, y, z, D, yon, F)
% 	float group_velocity(int x, int y, int z, float**** D, float* yon, float F)
% {
%   float detD;
%   float d11 = D[0][z][y][x];
%   float d12 = D[1][z][y][x];
%   float d13 = D[2][z][y][x];
%   float d22 = D[3][z][y][x];
%   float d23 = D[4][z][y][x];
%   float d33 = D[5][z][y][x];
%   detD = d11*(d22*d33 - pow(d23,2)) - d12*(d12*d33 - d13*d23) + d13*(d12*d23 - d22*d13);
%   //cout<<"determinant:"<<detD<<endl;

%   float x_tilda[3];
%   float X[3];
%   x_tilda[0] = yon[0]*(d22*d33-pow(d23,2))/detD + yon[1]*(d13*d23-d12*d33)/detD + yon[2]*(d12*d23 - d22*d13)/detD;
%   x_tilda[1] = yon[0]*(d23*d13-d12*d33)/detD + yon[1]*(d11*d33-pow(d13,2))/detD + yon[2]*(d13*d12 - d11*d23)/detD;
%   x_tilda[2] = yon[0]*(d12*d23-d22*d13)/detD + yon[1]*(d12*d13-d11*d23)/detD + yon[2]*(d11*d22 - pow(d12,2))/detD;
%   //cout<<"(x,y,z):"<<"("<<x_tilda[0]<<","<<x_tilda[1]<<","<<x_tilda[2]<<")"<<endl;
%   float A = sqrt( d11*pow(x_tilda[0],2) + d22*pow(x_tilda[1],2) + d33*pow(x_tilda[2],2) + 2*d12*x_tilda[0]*x_tilda[1] + 2*d13*x_tilda[0]*x_tilda[2] + 2*d23*x_tilda[1]*x_tilda[2])*F;
%   X[0] = x_tilda[0]/A;
%   X[1] = x_tilda[1]/A;
%   X[2] = x_tilda[2]/A;
%   float B = sqrt(d11*pow(X[0],2) + d22*pow(X[1],2) + d33*pow(X[2],2) + 2*d12*X[0]*X[1] + 2*d13*X[0]*X[2] + 2*d23*X[1]*X[2]);
  
%   return(sqrt(pow(X[0]*d11 + X[1]*d12 + X[2]*d13,2) + pow(X[0]*d12 + X[1]*d22 + X[2]*d23,2) + pow(X[0]*d13 + X[1]*d23 + X[2]*d33,2))*F/B);

% }
	d11 = D(1,z,y,x);
	d12 = D(2,z,y,x);
	d13 = D(3,z,y,x);
	d22 = D(4,z,y,x);
	d23 = D(5,z,y,x);
	d33 = D(6,z,y,x);
	vecI = [d11, d12, d13, d22, d23, d33];
	hessianmat = hessianvaluetomat(vecI);
	
	% The following code is the original implementation of determinant
	% detD = d11 * (d22 * d33 - d23 * d23) - d12 * (d12 * d33 - d13 * d23) + d13 * (d12 * d23 - d22 * d13);
	% fprintf('the determinant calculated by afm is %d\n', detD);
	detDmat = det(hessianmat);
	
	% fprintf('the determinant calculated by matlab is %d\n', detDmat);
	% x_tilda(1) = yon(1)*(d22 * d33 - d23 * d23)  / detDmat + yon(2) * (d13 * d23 - d12 * d33) / detDmat + yon(3) * (d12 * d23 - d22 * d13) / detDmat;
	% x_tilda(2) = yon(1)*(d23 * d13 - d12 * d33)  / detDmat + yon(2) * (d11 * d33 - d13 * d13) / detDmat + yon(3) * (d13 * d12 - d11 * d23) / detDmat;
	% x_tilda(3) = yon(1)*(d12 * d23 - d22 * d13)  / detDmat + yon(2) * (d12 * d13 - d11 * d23) / detDmat + yon(3) * (d11 * d22 - d12 * d12) / detDmat;
	% fprintf('tilda value 1: %d tilda value 2: %d tilda value 3: %d \n', x_tilda(1), x_tilda(2), x_tilda(3));
	
	% inv(hessianmat) : 3 x 3 yon : 3 x 1
	x_tildamat = inv(hessianmat) * yon;
	% fprintf('tilda value 1 matlab: %f tilda value 2 matlab : %f tilda value 3matlab : %f \n', x_tildamat(1), x_tildamat(2), x_tildamat(3));

	% x_tildamat : 3 x 1 hessianmat 3 x 3 A = (1 x 3) (3 x 3) (3 x 1) (1)
	% x_tildamat = [1;4;9]; test case
	% A = sqrt(d11*x_tildamat(1)*x_tildamat(1) + d22*x_tildamat(2)*x_tildamat(2) + d33*x_tildamat(3)*x_tildamat(3) + 2*d12*x_tildamat(1)*x_tildamat(2) + 2*d13*x_tildamat(1)*x_tildamat(3) + 2*d23*x_tildamat(2)*x_tildamat(3))*F;
	% fprintf('the value of A calculated by afm is %f\n', A)
	checksqrt = x_tildamat' * hessianmat * x_tildamat;
	if (checksqrt < 0)
		fprintf('You can not put negative number under the square happened in groupvelocity\n');
		% return;
        % This decision might revised in the future
        checksqrt = abs(checksqrt);
	end
	Amat = sqrt(checksqrt) * F;
	% fprintf('the value of A calculated by matlab is %f\n', Amat);	
	% X(1) = x_tildamat(1) / Amat;
	% X(2) = x_tildamat(2) / Amat;
	% X(3) = x_tildamat(3) / Amat;
	% X
	Xmat = x_tildamat / Amat;
	% Xmat
	% B = sqrt(d11*Xmat(1)*Xmat(1) + d22*Xmat(2)*Xmat(2) + d33*Xmat(3)*Xmat(3) + 2*d12*Xmat(1)*Xmat(2) + 2*d13*Xmat(1)*Xmat(3) + 2*d23*Xmat(2)*Xmat(3));
	% B
	Bmat = sqrt(Xmat' * hessianmat * Xmat);
	% Bmat
	% groupvvalue = sqrt((Xmat(1)*d11 + Xmat(2)*d12 + Xmat(3)*d13)^2 + (Xmat(1)*d12 + Xmat(2)*d22 + Xmat(3)*d23)^2 + (Xmat(1)*d13 + Xmat(2)*d23 + Xmat(3)*d33)^2)*F/Bmat;
	% groupvvalue
	temp = (hessianmat * Xmat).^2;
	groupvvaluemat = sqrt(sum(temp(:))) * F / Bmat;
	% groupvvaluemat
end



