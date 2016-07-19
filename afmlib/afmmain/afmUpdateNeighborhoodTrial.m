function [Tvalue, Ttag, value_iter] = afmUpdateNeighborhoodTrial(Tvalue, Ttag, Boundary, dx, dy, dz, afmSize, D, x, y, z, trX, trY, trZ, tetNo)
% void UpdateNeighborhoodTrial( float*** Tvalue, 
% 		int*** Ttag, float*** F, 
% 		unsigned char*** Boundary, 
% 		float dx, float dy, float dz, 
% 		int* Size, float**** D, int x, 
% 		int y, int z, harita* fm_map, 
% 		float** trX, float** trY, 
% 		float** trZ, int tetNo, 
% 		multimap<float,long>::iterator**** value_iter)
% {

	fprintf('afmUpdateNeighborhoodTrial function is created\n');

% 	float temp = 5000; // used to find the minimum value in the trial list.
% 	float temp2 = 5000; // used to find the minimum value in the trialC list.
	temp = 5000;
	temp2 = 5000;
	for(trNo = 1 : 8)
		% trNo
		% initialize flags for the roots and interest points
		flag1 = false;
		flag2 = false;
		flaga = true;
		flagb = true;
		flagc = true;
		% setting vectors
		a = [-1*trX(trNo,1), -1*trY(trNo,1), -1*trZ(trNo,1)];
		b = [-1*trX(trNo,2), -1*trY(trNo,2), -1*trZ(trNo,2)];
		c = [-1*trX(trNo,3), -1*trY(trNo,3), -1*trZ(trNo,3)];
		% a = [1 3 4] test case
		% tic
		d_a = sqrt((a(1)*dx)^2 + (a(2)*dy)^2 + (a(3)*dz)^2);
		d_b = sqrt((b(1)*dx)^2 + (b(2)*dy)^2 + (b(3)*dz)^2);
		d_c = sqrt((c(1)*dx)^2 + (c(2)*dy)^2 + (c(3)*dz)^2);
		% toc
		% tic
		% d_amat = sqrt(sum((a.*[dx dy dz]).^2))
		% toc
		a(1) = a(1)*dx/d_a;
		a(2) = a(2)*dy/d_a;
		a(3) = a(3)*dz/d_a;

		b(1) = b(1)*dx/d_b;
		b(2) = b(2)*dy/d_b;
		b(3) = b(3)*dz/d_b;

		c(1) = c(1)*dx/d_c;
		c(2) = c(2)*dy/d_c;
		c(3) = c(3)*dz/d_c;
		[p, q, r] = afmMatrixInverse3by3(a, b, c);
		% x = 3 test case
		conditionone = (x+trX(trNo,1))<1 || (y+trY(trNo,1))<1 || (z+trZ(trNo,1))<1; 
		% x = 23 test case
		conditiontwo = x+trX(trNo,1)>=afmSize(1) || y+trY(trNo,1)>= afmSize(2) || z+trZ(trNo,1)>= afmSize(3);
		conditionthree = (Boundary((z+trZ(trNo,1)), (y+trY(trNo,1)), (x+trX(trNo,1))) == 1);
        (z+trZ(trNo,1));
		(y+trY(trNo,1));
		(x+trX(trNo,1));
		indexxtemp = (x+trX(trNo,1));
		indexytemp = (y+trY(trNo,1));
		indexztemp = (z+trZ(trNo,1));
		conditionfour = (Ttag((z+trZ(trNo,1)),(y+trY(trNo,1)),(x+trX(trNo,1))) ~= 200);

		finalcondition = conditionone || conditiontwo || conditionthree || conditionfour;  
		if  finalcondition 
			flaga = false;
			Ta = 5000;
		else
			Ta = Tvalue((z+trZ(trNo,1)), (y+trY(trNo,1)), (x+trX(trNo,1)));
		end

		conditiontwo = x+trX(trNo,2)>=afmSize(1) || y+trY(trNo,2)>= afmSize(2) || z+trZ(trNo,2)>= afmSize(3);
		conditionthree = (Boundary((z+trZ(trNo,2)), (y+trY(trNo,2)), (x+trX(trNo,2))) == 1);
		conditionfour = (Ttag((z+trZ(trNo,2)), (y+trY(trNo,2)), (x+trX(trNo,2))) ~= 200);
		finalcondition = conditionone || conditiontwo || conditionthree || conditionfour;  
		if  finalcondition 
			flagb = false;
			Tb = 5000;
		else
			Tb = Tvalue((z+trZ(trNo,1)), (y+trY(trNo,1)), (x+trX(trNo,1)));
		end

		conditiontwo = x+trX(trNo,3)>=afmSize(1) || y+trY(trNo,3)>= afmSize(2) || z+trZ(trNo,3)>= afmSize(3);
		conditionthree = (Boundary((z+trZ(trNo,3)), (y+trY(trNo,3)), (x+trX(trNo,3))) == 1);
		conditionfour = (Ttag((z+trZ(trNo,3)), (y+trY(trNo,3)), (x+trX(trNo,3))) ~= 200);
		finalcondition = conditionone || conditiontwo || conditionthree || conditionfour;  
		if  finalcondition 
			flagc = false;
			Tc = 5000;
		else
			Tc = Tvalue((z+trZ(trNo,3)), (y+trY(trNo,3)), (x+trX(trNo,3)));
		end

		%if both interest points are ok.
		flaga = true; flagb = true; flagc = true;
		if(flaga && flagb && flagc)
			% Kx = (p(1)*Ta/d_a + p(2)*Tb/d_b + p(3)*Tc/d_c)
			% Ky = (q(1)*Ta/d_a + q(2)*Tb/d_b + q(3)*Tc/d_c)
			% Kz = (r(1)*Ta/d_a + r(2)*Tb/d_b + r(3)*Tc/d_c)
			pqrmat = [p;q;r];
			Kvec = pqrmat * [Ta/d_a; Tb/d_b; Tc/d_c];
			Kx = Kvec(1);
			Ky = Kvec(2);
			Kz = Kvec(3);
			% Cx = (p(1)/d_a + p(2)/d_b + p(3)/d_c)
			% Cy = (q(1)/d_a + q(2)/d_b + q(3)/d_c)
			% Cz = (r(1)/d_a + r(2)/d_b + r(3)/d_c)
			Cmat = pqrmat * [1/d_a; 1/d_b; 1/d_c];
			Cx = Cmat(1);
			Cy = Cmat(2);
			Cz = Cmat(3);
			size(Cmat)
			w1 = D[0,z,y,x]*pow(Cx,2) + D[3,z,y,x]*pow(Cy,2) + D[5,z,y,x]*pow(Cz,2) + 2*D[1,z,y,x]*Cx*Cy + 2*D[2,z,y,x]*Cx*Cz + 2*D[4][z][y][x]*Cy*Cz;
			% Dhessianvec = [D(0,z,y,x), D(1,z,y,x), D(2,z,y,x), D(3,z,y,x), D(4,z,y,x), D(5,z,y,x)];
			% DHmat = hessianvaluetomat(Dhessianvec);

		end

	end
% 	int trNo;
% 	float a[3], b[3], c[3], Ta, Tb, Tc, p[3], q[3], 
% 				r[3], Cx, Cy, Cz, Kx, Ky, Kz, w1, w2, w3, R[2], 
% 				dT[3], d[3], dir1, dir2, dir3, vga, vgb, vgc, 
% 				d_a, d_b, d_c;
% 	bool flag1, flag2, flaga, flagb, flagc, flag_imag;

% 	multimap<float,long>::iterator viter;
% 	long index;
% 	//mexPrintf("mex print is working");
% 	for(trNo=0; trNo < tetNo; trNo++)
% 	{
% 		// initializing flags for the roots and interest points
% 		flag1 = false;
% 		flag2 = false;
% 		flaga = true;
% 		flagb = true;
% 		flagc = true;
% 		// setting vectors.
% 		a[0] = -1*trX[trNo][0];
% 		a[1] = -1*trY[trNo][0];
% 		a[2] = -1*trZ[trNo][0];

% 		b[0] = -1*trX[trNo][1];
% 		b[1] = -1*trY[trNo][1];
% 		b[2] = -1*trZ[trNo][1];

% 		c[0] = -1*trX[trNo][2];
% 		c[1] = -1*trY[trNo][2];
% 		c[2] = -1*trZ[trNo][2];

% 		d_a = sqrt(pow(a[0]*dx,2) + pow(a[1]*dy,2) + pow(a[2]*dz,2));
% 		d_b = sqrt(pow(b[0]*dx,2) + pow(b[1]*dy,2) + pow(b[2]*dz,2));
% 		d_c = sqrt(pow(c[0]*dx,2) + pow(c[1]*dy,2) + pow(c[2]*dz,2));
% 		// checking if interest points are ok.
% 		a[0] = a[0]*dx/d_a;
% 		a[1] = a[1]*dy/d_a;
% 		a[2] = a[2]*dz/d_a;
% 		b[0] = b[0]*dx/d_b;
% 		b[1] = b[1]*dy/d_b;
% 		b[2] = b[2]*dz/d_b;
% 		c[0] = c[0]*dx/d_c;
% 		c[1] = c[1]*dy/d_c;
% 		c[2] = c[2]*dz/d_c;
% 		MatrixInverse3by3( a,b,c,p,q,r );

% 		if(x+trX[trNo][0]<0 || y+trY[trNo][0]<0 || z+trZ[trNo][0]<0 || x+trX[trNo][0]>=Size[0] || y+trY[trNo][0]>=Size[1] || z+trZ[trNo][0]>=Size[2] || Boundary[(unsigned short)(z+trZ[trNo][0])][(unsigned short)(y+trY[trNo][0])][(unsigned short)(x+trX[trNo][0])] == 1 || Ttag[(unsigned short)(z+trZ[trNo][0])][(unsigned short)(y+trY[trNo][0])][(unsigned short)(x+trX[trNo][0])] != 200)
% 		{
% 			flaga = false;
% 			Ta = 5000;
% 		}
% 		else
% 		{
% 			Ta = Tvalue[(unsigned short)(z+trZ[trNo][0])][(unsigned short)(y+trY[trNo][0])][(unsigned short)(x+trX[trNo][0])];
% 		}

% 		if(x+trX[trNo][1]<0 || y+trY[trNo][1]<0 || z+trZ[trNo][1]<0 || x+trX[trNo][1]>=Size[0] || y+trY[trNo][1]>=Size[1] || z+trZ[trNo][1]>=Size[2] || Boundary[(unsigned short)(z+trZ[trNo][1])][(unsigned short)(y+trY[trNo][1])][(unsigned short)(x+trX[trNo][1])] == 1 || Ttag[(unsigned short)(z+trZ[trNo][1])][(unsigned short)(y+trY[trNo][1])][(unsigned short)(x+trX[trNo][1])] != 200)
% 		{
% 			flagb = false;
% 			Tb = 5000;
% 		}
% 		else
% 		{
% 			Tb = Tvalue[(unsigned short)(z+trZ[trNo][1])][(unsigned short)(y+trY[trNo][1])][(unsigned short)(x+trX[trNo][1])];
% 		}

% 		if(x+trX[trNo][2]<0 || y+trY[trNo][2]<0 || z+trZ[trNo][2]<0 || x+trX[trNo][2]>=Size[0] || y+trY[trNo][2]>=Size[1] || z+trZ[trNo][2]>=Size[2] || Boundary[(unsigned short)(z+trZ[trNo][2])][(unsigned short)(y+trY[trNo][2])][(unsigned short)(x+trX[trNo][2])] == 1 || Ttag[(unsigned short)(z+trZ[trNo][2])][(unsigned short)(y+trY[trNo][2])][(unsigned short)(x+trX[trNo][2])] != 200)
% 		{
% 			flagc = false;
% 			Tc = 5000;
% 		}
% 		else
% 		{
% 			Tc = Tvalue[(unsigned short)(z+trZ[trNo][2])][(unsigned short)(y+trY[trNo][2])][(unsigned short)(x+trX[trNo][2])];
% 		}



% 		// if both interest points are ok.
% 		if(flaga && flagb && flagc)
% 		{
% 			Kx = (p[0]*Ta/d_a + p[1]*Tb/d_b + p[2]*Tc/d_c);
% 			Ky = (q[0]*Ta/d_a + q[1]*Tb/d_b + q[2]*Tc/d_c);
% 			Kz = (r[0]*Ta/d_a + r[1]*Tb/d_b + r[2]*Tc/d_c);
% 			Cx = (p[0]/d_a + p[1]/d_b + p[2]/d_c);
% 			Cy = (q[0]/d_a + q[1]/d_b + q[2]/d_c);
% 			Cz = (r[0]/d_a + r[1]/d_b + r[2]/d_c);

% 			w1 = D[0][z][y][x]*pow(Cx,2) + D[3][z][y][x]*pow(Cy,2) + D[5][z][y][x]*pow(Cz,2) + 2*D[1][z][y][x]*Cx*Cy + 2*D[2][z][y][x]*Cx*Cz + 2*D[4][z][y][x]*Cy*Cz;


% 			w2 = -2*D[0][z][y][x]*Cx*Kx - 2*D[3][z][y][x]*Cy*Ky - 2*D[5][z][y][x]*Cz*Kz - 2*D[1][z][y][x]*Cx*Ky - 2*D[1][z][y][x]*Cy*Kx - 2*D[2][z][y][x]*Cx*Kz - 2*D[2][z][y][x]*Cz*Kx - 2*D[4][z][y][x]*Cy*Kz - 2*D[4][z][y][x]*Cz*Ky;

% 			w3 = D[0][z][y][x]*pow(Kx,2) + D[3][z][y][x]*pow(Ky,2) + D[5][z][y][x]*pow(Kz,2) + 2*D[1][z][y][x]*Kx*Ky + 2*D[2][z][y][x]*Kx*Kz + 2*D[4][z][y][x]*Ky*Kz - 1/pow(F[z][y][x],2);


% 			flag_imag = find_roots_indic( R,100*w1,100*w2,100*w3 );
% 			// doing the update regularly if the roots are ok. 
% 			if(!flag_imag && R[0] > 0)
% 			{
% 				dT[0] = p[0]*(R[0] - Ta)/d_a + p[1]*(R[0] - Tb)/d_b + p[2]*(R[0] - Tc)/d_c;
% 				dT[1] = q[0]*(R[0] - Ta)/d_a + q[1]*(R[0] - Tb)/d_b + q[2]*(R[0] - Tc)/d_c;
% 				dT[2] = r[0]*(R[0] - Ta)/d_a + r[1]*(R[0] - Tb)/d_b + r[2]*(R[0] - Tc)/d_c;

% 				d[0] = (D[0][z][y][x]*dT[0] + D[1][z][y][x]*dT[1] + D[2][z][y][x]*dT[2]);
% 				d[1] = (D[1][z][y][x]*dT[0] + D[3][z][y][x]*dT[1] + D[4][z][y][x]*dT[2]);
% 				d[2] = (D[2][z][y][x]*dT[0] + D[4][z][y][x]*dT[1] + D[5][z][y][x]*dT[2]);

% 				dir1 = (d[0]*(a[1]*b[2] - a[2]*b[1]) - d[1]*(a[0]*b[2] - a[2]*b[0]) + d[2]*(a[0]*b[1] - a[1]*b[0]));
% 				dir2 = (d[0]*(b[1]*c[2] - b[2]*c[1]) - d[1]*(b[0]*c[2] - b[2]*c[0]) + d[2]*(b[0]*c[1] - b[1]*c[0]));
% 				dir3 = (d[0]*(c[1]*a[2] - c[2]*a[1]) - d[1]*(c[0]*a[2] - c[2]*a[0]) + d[2]*(c[0]*a[1] - c[1]*a[0]));

% 				if(( (sign(dir1)<0 && sign(dir2)<=0 && sign(dir3)<=0) || (sign(dir1)<=0 && sign(dir2)<=0 && sign(dir3)<0) || (sign(dir1)<=0 && sign(dir2)<0 && sign(dir3)<=0))&& 1)
% 				{
% 					temp = min( temp,R[0] );
% 					flag1 = true;
% 				}
% 			}

% 			if(!flag_imag && R[1] > 0)
% 			{
% 				dT[0] = p[0]*(R[1] - Ta)/d_a + p[1]*(R[1] - Tb)/d_b + p[2]*(R[1] - Tc)/d_c;
% 				dT[1] = q[0]*(R[1] - Ta)/d_a + q[1]*(R[1] - Tb)/d_b + q[2]*(R[1] - Tc)/d_c;
% 				dT[2] = r[0]*(R[1] - Ta)/d_a + r[1]*(R[1] - Tb)/d_b + r[2]*(R[1] - Tc)/d_c;

% 				d[0] = (D[0][z][y][x]*dT[0] + D[1][z][y][x]*dT[1] + D[2][z][y][x]*dT[2]);
% 				d[1] = (D[1][z][y][x]*dT[0] + D[3][z][y][x]*dT[1] + D[4][z][y][x]*dT[2]);
% 				d[2] = (D[2][z][y][x]*dT[0] + D[4][z][y][x]*dT[1] + D[5][z][y][x]*dT[2]);

% 				dir1 = (d[0]*(a[1]*b[2] - a[2]*b[1]) - d[1]*(a[0]*b[2] - a[2]*b[0]) + d[2]*(a[0]*b[1] - a[1]*b[0]));
% 				dir2 = (d[0]*(b[1]*c[2] - b[2]*c[1]) - d[1]*(b[0]*c[2] - b[2]*c[0]) + d[2]*(b[0]*c[1] - b[1]*c[0]));
% 				dir3 = (d[0]*(c[1]*a[2] - c[2]*a[1]) - d[1]*(c[0]*a[2] - c[2]*a[0]) + d[2]*(c[0]*a[1] - c[1]*a[0]));

% 				if(( (sign(dir1)<0 && sign(dir2)<=0 && sign(dir3)<=0) 
% 							|| (sign(dir1)<=0 && sign(dir2)<=0 && sign(dir3)<0) 
% 							|| (sign(dir1)<=0 && sign(dir2)<0 && sign(dir3)<=0))&& 1)
% 				{
% 					temp = min( temp,R[1] );
% 					flag2 = true;
% 				}
% 			}
% 		}

% 		// if found roots are not useful. both flags should be false.
% 		if( (!flag1 && !flag2) || 1)
% 		{

% 			if(flaga && flagb && flagc) // all interest points are ok.
% 			{

% 				temp = min(temp, minimize_Analytic2(Ta,Tb,x,y,z,D,a,b,F[z][y][x],d_a,d_b)  );

% 				temp = min(temp, minimize_Analytic2(Ta,Tc,x,y,z,D,a,c,F[z][y][x],d_a,d_c)  );

% 				temp = min(temp, minimize_Analytic2(Tb,Tc,x,y,z,D,b,c,F[z][y][x],d_b,d_c)  );

% 			}	
% 			else if(!flaga && flagb && flagc) // only two are ok.
% 			{

% 				temp = min(temp, minimize_Analytic2(Tb,Tc,x,y,z,D,b,c,F[z][y][x],d_b,d_c)  );

% 			}
% 			else if(flaga && !flagb && flagc) // only two are ok.
% 			{

% 				temp = min(temp, minimize_Analytic2(Ta,Tc,x,y,z,D,a,c,F[z][y][x],d_a,d_c)  );

% 			}
% 			else if(flaga && flagb && !flagc) // only two are ok.
% 			{

% 				temp = min(temp, minimize_Analytic2(Ta,Tb,x,y,z,D,a,b,F[z][y][x],d_a,d_b)  );

% 			}
% 			else if(flaga && !flagb && !flagc) // only one is ok.
% 			{

% 				vga = group_velocity(x,y,z,D,a,F[z][y][x]);

% 				temp = min( temp, Ta+d_a/vga );

% 			}
% 			else if(!flaga && flagb && !flagc) // only one is ok.
% 			{
% 				vgb = group_velocity(x,y,z,D,b,F[z][y][x]);

% 				temp = min( temp, Tb+d_b/vgb );

% 			}
% 			else if(!flaga && !flagb && flagc) // only one is ok.
% 			{			    			 
% 				vgc = group_velocity(x,y,z,D,c,F[z][y][x]);	

% 				temp = min( temp, Tc+d_c/vgc );

% 			}
% 		}

% 	}


% 	if(temp != 5000)
% 	{


% 		if( Ttag[z][y][x] == 175 ) // point is already in the good trial list. UPDATE.
% 		{

% 			viter = value_iter[z][y][x][0];
% 			if( viter->first > temp )
% 			{
% 				Tvalue[z][y][x] = temp;
% 				index = sub2ind(x,y,z,Size);
% 				fm_map->trial.erase( viter );
% 				value_iter[z][y][x][0] = fm_map->trial.insert( make_pair(temp,index) );
% 			}
% 		}

% 		else if( Ttag[z][y][x] == 125 ) // point is not in the good list but is in the bad list. REMOVE AND ADD TO THE GOOD LIST.
% 		{

% 			if( Tvalue[z][y][x] >= temp )
% 			{

% 				// removing from the bad one.
% 				viter = value_iter[z][y][x][1];	
% 				fm_map->trialC.erase( viter );

% 				// adding to the good trial list
% 				index = sub2ind(x,y,z,Size);
% 				value_iter[z][y][x][0] = fm_map->trial.insert( make_pair(temp,index) ); //  write the address.
% 				Tvalue[z][y][x] = temp;
% 				Ttag[z][y][x] = 175; // making good trial tag.

% 			}

% 		}
% 		else // it is not in both lists
% 		{	

% 			index = sub2ind(x,y,z,Size);
% 			value_iter[z][y][x][0] = fm_map->trial.insert( make_pair(temp,index) ); // write the address.
% 			Tvalue[z][y][x] = temp;
% 			Ttag[z][y][x] = 175; // good list trial number

% 		}



% 	}

% 	if(temp2 != 5000)
% 	{
% 		// distance to trial lists.

% 		if( Ttag[z][y][x] == 175 ) // point is alread in the good trial list. CHECK-REMOVE-UDPATE.
% 		{

% 			if(Tvalue[z][y][x] > temp2) // found a better candidate.
% 			{

% 				// removing from the good trial list.
% 				viter = value_iter[z][y][x][0];
% 				fm_map->trial.erase(viter);

% 				// adding to the bad trial list.
% 				index = sub2ind(x,y,z,Size);
% 				value_iter[z][y][x][1] = fm_map->trialC.insert( make_pair(temp2,index) ); // write the address
% 				Tvalue[z][y][x] = temp2;
% 				Ttag[z][y][x] = 125; // making bad trial tag.

% 			}
% 		}		

% 		else if( Ttag[z][y][x] == 125 ) // not in the good list but in the bad list.
% 		{

% 			viter = value_iter[z][y][x][1];
% 			if(viter->first > temp2)
% 			{
% 				Tvalue[z][y][x] = temp2;
% 				index = sub2ind(x,y,z,Size);
% 				fm_map->trialC.erase(viter);
% 				value_iter[z][y][x][1] = fm_map->trialC.insert( make_pair(temp2,index) );
% 			}

% 		}
% 		else // it is not in both lists
% 		{

% 			index = sub2ind(x,y,z,Size);
% 			value_iter[z][y][x][1] = fm_map->trialC.insert( make_pair(temp2,index) ); // write the address
% 			Tvalue[z][y][x] = temp2;
% 			Ttag[z][y][x] = 125; //  bad list trial number


% 		}

% 	}


% }
	Tvalue = [];
	Ttag = [];
	value_iter = [];
end
