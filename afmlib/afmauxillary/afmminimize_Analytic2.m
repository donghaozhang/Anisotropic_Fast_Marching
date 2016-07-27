function f = minimize_Analytic2(Ta, Tb, x, y, z, Dbig, a, b, F, d_a, d_b)
% float minimize_Analytic2( float Ta, float Tb, int x, int y, int z
% 			  , float**** Dbig, float* a, float* b, float F
% 			  ,float d_a, float d_b )
% {
  
%   bool flag_imag;
%   float t, A1, A2, B1, B2, C1, C2, f_atZero, f_atOne, f_q1, f_q2, f, detD;
%   float temp1,temp2;
%   float w[3], q[2], P[3][3], vec[3];
  
  
%   float d11 = Dbig[0][z][y][x];
%   float d12 = Dbig[1][z][y][x];
%   float d13 = Dbig[2][z][y][x];
%   float d22 = Dbig[3][z][y][x];
%   float d23 = Dbig[4][z][y][x];
%   float d33 = Dbig[5][z][y][x];
  
  d11 = Dbig(1,z,y,x);
  d12 = Dbig(2,z,y,x);
  d13 = Dbig(3,z,y,x);
  d22 = Dbig(4,z,y,x);
  d23 = Dbig(5,z,y,x);
  d33 = Dbig(6,z,y,x);
  vecI = [d11, d12, d13, d22, d23, d33];
%   detD = d11*(d22*d33 - std::pow(d23,2)) - d12*(d12*d33 - d13*d23) + 
%     d13*(d12*d23 - d22*d13);
  P = hessianvaluetomat(vecI);
  detDmat = det(P);  
  % detD = d11 * (d22 * d33 - d23 * d23) - d12 * (d12 * d33 - d13 * d23) + d13 * (d12 * d23 - d22 * d13);  

%   P[0][0] = (d22*d33-std::pow(d23,2))/detD;
%   P[0][1] = (d23*d13-d12*d33)/detD;
%   P[0][2] = (d12*d23-d22*d13)/detD;
%   P[1][0] = (d13*d23-d12*d33)/detD;
%   P[1][1] = (d11*d33-std::pow(d13,2))/detD;
%   P[1][2] = (d12*d13-d11*d23)/detD;
%   P[2][0] = (d12*d23 - d22*d13)/detD;
%   P[2][1] = (d13*d12 - d11*d23)/detD;
%   P[2][2] = (d11*d22 - std::pow(d12,2))/detD;

  % P(1, 1) = (d22 * d33 - d23 * d23) / detDmat;
  % P(1, 2) = (d23 * d13 - d12 * d33) / detDmat;
  % P(1, 3) = (d12 * d23 - d22 * d13) / detDmat;
  % P(2, 1) = (d13 * d23 - d12 * d33) / detDmat;
  % P(2, 2) = (d11 * d33 - d13 * d13) / detDmat;
  % P(2, 3) = (d12 * d13 - d11 * d23) / detDmat;
  % P(3, 1) = (d12 * d23 - d22 * d13) / detDmat;
  % P(3, 2) = (d13 * d12 - d11 * d23) / detDmat;
  % P(3, 3) = (d11 * d22 - d12 * d12) / detDmat;
  % P
  matP = inv(P);
  % matP
  
%   w[0] = d_a*a[0] - d_b*b[0];
%   w[1] = d_a*a[1] - d_b*b[1];
%   w[2] = d_a*a[2] - d_b*b[2];
  % w(1) = d_a * a(1) - d_b * b(1);
  % w(2) = d_a * a(2) - d_b * b(2);
  % w(3) = d_a * a(3) - d_b * b(3);
  % w
  wmat = d_a * a - d_b * b;
  % wmat
  % size(wmat)
  % wmat
  
%   A1 = (P[0][0]*std::pow(a[0],2) + P[1][1]*std::pow(a[1],2) + P[2][2]*std::pow(a[2],2))*std::pow(d_a,2) + (P[0][0]*std::pow(b[0],2) + P[1][1]*std::pow(b[1],2) + P[2][2]*std::pow(b[2],2))*std::pow(d_b,2) - 2*d_a*d_b*(P[0][0]*a[0]*b[0] + P[1][1]*a[1]*b[1] + P[2][2]*a[2]*b[2]) + 2*std::pow(d_a,2)*(P[0][1]*a[0]*a[1] + P[0][2]*a[0]*a[2] + P[1][2]*a[1]*a[2]) + 2*std::pow(d_b,2)*(P[0][1]*b[0]*b[1] + P[0][2]*b[0]*b[2] + P[1][2]*b[1]*b[2]) - 2*d_a*d_b*(P[0][1]*a[0]*b[1] + P[0][2]*a[0]*b[2] + P[1][2]*a[1]*b[2]) - 2*d_a*d_b*(P[0][1]*b[0]*a[1] + P[0][2]*b[0]*a[2] + P[1][2]*b[1]*a[2]);
   % A1 = (P(1,1)*(a(1)^2) + P(2,2)*(a(2)^2) + P(3,3)*(a(3)^2))*(d_a^2) + (P(1,1)*(b(1)^2) + P(2,2)*(b(2)^2) + P(3,3)*(b(3)^2))*(d_b^2) - 2*d_a*d_b*(P(1,1)*a(1)*b(1) + P(2,2)*a(2)*b(2) + P(3,3)*a(3)*b(3)) + 2*(d_a^2)*(P(1,2)*a(1)*a(2) + P(1,3)*a(1)*a(3) + P(2,3)*a(2)*a(3)) + 2*(d_b^2)*(P(1,2)*b(1)*b(2) + P(1,3)*b(1)*b(3) + P(2,3)*b(2)*b(3)) - 2*d_a*d_b*(P(1,2)*a(1)*b(2) + P(1,3)*a(1)*b(3) + P(2,3)*a(2)*b(3)) - 2*d_a*d_b*(P(1,2)*b(1)*a(2) + P(1,3)*b(1)*a(3) + P(2,3)*b(2)*a(3))
   A1mat = d_a^2 * a * P * a' + d_b^2 * b * P * b'- 2 * d_a * d_b * a * P * b';
%   B1 = 2*d_a*d_b*(P[0][0]*a[0]*b[0] + P[1][1]*a[1]*b[1] + P[2][2]*a[2]*b[2]) - 2*std::pow(d_b,2)*(P[0][0]*std::pow(b[0],2) + P[1][1]*std::pow(b[1],2) + P[2][2]*std::pow(b[2],2)) + 2*d_a*d_b*(P[0][1]*a[0]*b[1] + P[0][2]*a[0]*b[2] + P[1][2]*a[1]*b[2]) + 2*d_a*d_b*(P[0][1]*b[0]*a[1] + P[0][2]*b[0]*a[2] + P[1][2]*b[1]*a[2]) - 4*std::pow(d_b,2)*(P[0][1]*b[0]*b[1] + P[0][2]*b[0]*b[2] + P[1][2]*b[1]*b[2]);
  % B1 = 2*d_a*d_b*(P(1,1)*a(1)*b(1) + P(2,2)*a(2)*b(2) + P(3,3)*a(3)*b(3)) - 2*(d_b^2)*(P(1,1)*(b(1)^2) + P(2,2)*(b(2)^2) + P(3,3)*(b(3)^2)) + 2*d_a*d_b*(P(1,2)*a(1)*b(2) + P(1,3)*a(1)*b(3) + P(2,3)*a(2)*b(3)) + 2*d_a*d_b*(P(1,2)*b(1)*a(2) + P(1,3)*b(1)*a(3) + P(2,3)*b(2)*a(3)) - 4*(d_b^2)*(P(1,2)*b(1)*b(2) + P(1,3)*b(1)*b(3) + P(2,3)*b(2)*b(3));
  % B1temp = 2*d_a*d_b*(P(1,1)*a(1)*b(1) + P(2,2)*a(2)*b(2) + P(3,3)*a(3)*b(3)) + 2*d_a*d_b*(P(1,2)*a(1)*b(2) + P(1,3)*a(1)*b(3) + P(2,3)*a(2)*b(3)) + + 2*d_a*d_b*(P(1,2)*b(1)*a(2) + P(1,3)*b(1)*a(3) + P(2,3)*b(2)*a(3))
  % B1
  B1mat = 2 * d_a * a * P * d_b * b' - 2 * d_b * d_b * b * P * b';
  % B1mat
  % B1mattemp =   2 * d_a * a * P * d_b * b' 

%   C1 = std::pow(d_b,2)*(P[0][0]*std::pow(b[0],2) + P[1][1]*std::pow(b[1],2) + P[2][2]*std::pow(b[2],2)) + 2*std::pow(d_b,2)*(P[0][1]*b[0]*b[1] + P[0][2]*b[0]*b[2] + P[1][2]*b[1]*b[2]);
  % b Pb b' = 1 x 3 3 x 3 3 x 1 
  % C1 = (d_b^2)*(P(1,1)*(b(1)^2) + P(2,2)*(b(2)^2) + P(3,3)*(b(3)^2)) + 2*(d_b^2)*(P(1,2)*b(1)*b(2) + P(1,3)*b(1)*b(3) + P(2,3)*b(2)*b(3));
  % C1
  C1mat = (d_b^2) * (b * P * b');
  % C1mat

%   temp1 = P[0][0]*std::pow(w[0],2) + P[1][1]*std::pow(w[1],2) + P[2][2]*std::pow(w[2],2) + 2*P[0][1]*w[0]*w[1] + 2*P[0][2]*w[0]*w[2] + 2*P[1][2]*w[1]*w[2];
  % temp1 = P(1,1)*(wmat(1)^2) + P(2,2)*(wmat(2)^2) + P(3,3)*(wmat(3)^2) + 2*P(1,2)*wmat(1)*wmat(2) + 2*P(1,3)*wmat(1)*wmat(3) + 2*P(2,3)*wmat(2)*wmat(3)
  temp1mat = wmat * P * wmat'; 
%   temp2 = d_b*( P[0][0]*w[0]*b[0] + P[1][1]*w[1]*b[1] + P[2][2]*w[2]*b[2] + P[0][1]*w[0]*b[1] + P[0][1]*w[1]*b[0] + P[0][2]*w[0]*b[2] + P[0][2]*w[2]*b[0] + P[1][2]*w[1]*b[2] + P[1][2]*w[2]*b[1] );
  % temp2 = d_b*( P(1,1)*wmat(1)*b(1) + P(2,2)*wmat(2)*b(2) + P(3,3)*wmat(3)*b(3) + P(1,2)*wmat(1)*b(2) + P(1,2)*wmat(2)*b(1) + P(1,3)*wmat(1)*b(3) + P(1,3)*wmat(3)*b(1) + P(2,3)*wmat(2)*b(3) + P(2,3)*wmat(3)*b(2) )
  temp2mat = wmat * P * b' * d_b;
%   A2 = std::pow(temp1,2);
  A2 = temp1mat^2;
  B2 = 2*temp1mat*temp2mat;
  % B2 = 2 * temp1 * temp2;
%   C2 = std::pow(temp2,2);
  C2 = temp2mat^2;
  
%   A1 = A1*std::pow(Tb - Ta,2)*std::pow(F,2);
  A1mat = A1mat * ((Tb - Ta)^2) * (F^2);
%   B1 = B1*std::pow(Tb - Ta,2)*std::pow(F,2);
  B1mat = B1mat * ((Tb - Ta)^2) * (F^2);  
%   C1 = C1*std::pow(Tb - Ta,2)*std::pow(F,2);
  C1mat = C1mat * ((Tb - Ta)^2) * (F^2);

  
  
  
  
%   flag_imag = find_roots_indic( q, 100*(A2-A1), 100*(B2-B1), 100*(C2-C1) );
  [q, deltalogic] = afmfind_roots_indic(100 * (A2 - A1mat), 100 * (B2 - B1mat), 100 * (C2 - C1mat));

%   if(flag_imag)
%   {
%   	q[0] = -1;
%   	q[1] = -1;
%   }
  if(deltalogic)    
      q(1) = -1;
      q(2) = -1;
  end
    
%   t = 0;
  t = 0;
%   vec[0] = a[0]*t*d_a + b[0]*(1-t)*d_b;
%   vec[1] = a[1]*t*d_a + b[1]*(1-t)*d_b;
%   vec[2] = a[2]*t*d_a + b[2]*(1-t)*d_b;
  % vec(1) = a(1) * t * d_a + b(1) * (1 - t) * d_b;
  % vec(2) = a(2) * t * d_a + b(2) * (1 - t) * d_b;
  % vec(3) = a(3) * t * d_a + b(3) * (1 - t) * d_b;
  % vec
  vecmat = a * t * d_a + b * (1 - t) * d_b;      
%   f_atZero = Ta*t + (1-t)*Tb + sqrt( P[0][0]*std::pow(vec[0],2) + P[1][1]*std::pow(vec[1],2) + P[2][2]*std::pow(vec[2],2) + 2*P[0][1]*vec[0]*vec[1] + 2*P[0][2]*vec[0]*vec[2] + 2*P[1][2]*vec[1]*vec[2] )/F;  
  % f_atZero = Ta*t + (1-t)*Tb + sqrt( P(1,1)*(vecmat(1)^2) + P(2,2)*(vecmat(2)^2) + P(3,3)*(vecmat(3)^2) + 2*P(1,2)*vecmat(1)*vecmat(2) + 2*P(1,3)*vecmat(1)*vecmat(3) + 2*P(2,3)*vecmat(2)*vecmat(3) )/F  
  f_atZeromat = Ta*t + (1-t)*Tb + sqrt(vecmat * P * vecmat') / F;
%   t = 1;
%   vec[0] = a[0]*t*d_a + b[0]*(1-t)*d_b;
%   vec[1] = a[1]*t*d_a + b[1]*(1-t)*d_b;
%   vec[2] = a[2]*t*d_a + b[2]*(1-t)*d_b;
  t = 1;
  vecmat = a * t * d_a + b * (1 - t) * d_b; 
  f_atOnemat = Ta*t + (1-t)*Tb + sqrt(vecmat * P * vecmat') / F;    
%   f_atOne = Ta*t + (1-t)*Tb + sqrt( P[0][0]*std::pow(vec[0],2) + P[1][1]*std::pow(vec[1],2) + P[2][2]*std::pow(vec[2],2) + 2*P[0][1]*vec[0]*vec[1] + 2*P[0][2]*vec[0]*vec[2] + 2*P[1][2]*vec[1]*vec[2] )/F;
%   f = min(f_atOne, f_atZero);
  f = afmmin(f_atZeromat, f_atOnemat);
  f_realflag = isreal(f);
  if ~f_realflag
    fprintf('Imaginary value appears stage one\n');
  end
%   if(q[0] >= 0 && q[0] <= 1)
%   {
%   	t = q[0];
%   	vec[0] = a[0]*t*d_a + b[0]*(1-t)*d_b;
%   	vec[1] = a[1]*t*d_a + b[1]*(1-t)*d_b;
%   	vec[2] = a[2]*t*d_a + b[2]*(1-t)*d_b;
      
%   	f_q1 = Ta*t + (1-t)*Tb + sqrt( P[0][0]*std::pow(vec[0],2) + P[1][1]*std::pow(vec[1],2) + P[2][2]*std::pow(vec[2],2) + 2*P[0][1]*vec[0]*vec[1] + 2*P[0][2]*vec[0]*vec[2] + 2*P[1][2]*vec[1]*vec[2] )/F;
%   	f = min(f, f_q1);
%   }
  % q(1) = 0.5; Just for testing purpose
  if(q(1) >= 0 && q(1) <= 1)
    t = q(1);
    % vec(1) = a(1)*t*d_a + b(1)*(1-t)*d_b;
    % vec(2) = a(2)*t*d_a + b(2)*(1-t)*d_b;
    % vec(3) = a(3)*t*d_a + b(3)*(1-t)*d_b;
    % vec
    vecmat = a * t * d_a + b * (1 - t) * d_b;
    % f_q1 = Ta*t + (1-t)*Tb + sqrt( P(1, 1)*(vecmat(1)^2) + P(2, 2)*(vecmat(2)^2) + P(3, 3)*(vecmat(3)^2) + 2*P(1, 2)*vecmat(1)*vecmat(2) + 2*P(1, 3)*vecmat(1)*vecmat(3) + 2*P(2, 3)*vecmat(2)*vecmat(3) )/F;
    f_q1mat = Ta*t + (1-t)*Tb + sqrt(vecmat * P * vecmat') / F;
    f = afmmin(f, f_q1mat);
    f_realflag = isreal(f);
    if ~f_realflag
      fprintf('Imaginary value appears stage two\n');
    end
  end
%   if(q[1] >= 0 && q[1] <= 1)
%   {
%   	t = q[1];
%   	vec[0] = a[0]*t*d_a + b[0]*(1-t)*d_b;
%   	vec[1] = a[1]*t*d_a + b[1]*(1-t)*d_b;
%   	vec[2] = a[2]*t*d_a + b[2]*(1-t)*d_b;
      
%   	f_q2 = Ta*t + (1-t)*Tb + sqrt( P[0][0]*std::pow(vec[0],2) + P[1][1]*std::pow(vec[1],2) + P[2][2]*std::pow(vec[2],2) + 2*P[0][1]*vec[0]*vec[1] + 2*P[0][2]*vec[0]*vec[2] + 2*P[1][2]*vec[1]*vec[2] )/F;
%   	f = min(f, f_q2);
%   }
  % q(2) = 0.5;
  if (q(2) >= 0 && q(2) <= 1)
    t = q(2);
    % vec(1) = a(1)*t*d_a + b(1)*(1-t)*d_b;
    % vec(2) = a(2)*t*d_a + b(2)*(1-t)*d_b;
    % vec(3) = a(3)*t*d_a + b(3)*(1-t)*d_b;
    % vec
    vecmat = a * t * d_a + b * (1 - t) * d_b;
    % f_q2 = Ta*t + (1-t)*Tb + sqrt( P(1,1)*(vecmat(1)^2) + P(2,2)*(vecmat(2)^2) + P(3,3)*(vecmat(3)^2) + 2*P(1,2)*vecmat(1)*vecmat(2) + 2*P(1,3)*vecmat(1)*vecmat(3) + 2*P(2,3)*vecmat(2)*vecmat(3) )/F
    f_q2mat = Ta*t + (1-t)*Tb + sqrt(vecmat * P * vecmat') / F;
    f = afmmin(f, f_q2mat);
  end
%   return(f);
    f_realflag = isreal(f);
    if ~f_realflag
      fprintf('Imaginary value appears stage two\n');
    end
    f
% }
end