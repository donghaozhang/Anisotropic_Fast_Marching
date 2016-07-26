function [Tvalue, Ttag, value_iter] = afmUpdateNeighborhoodTrial(Tvalue, Ttag, F, Boundary, dx, dy, dz, afmSize, D, x, y, z, trial, trialC, trX, trY, trZ, tetNo, value_iter)
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
		% !!!!!!! delete the following line in the future
		% flaga = true; flagb = true; flagc = true; %!!!!!!!!! delete it in the future
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
			% w1 = D(1,z,y,x)*(Cx^2) + D(4,z,y,x)*(Cy^2) + D(6,z,y,x)*(Cz^2) + 2*D(2,z,y,x)*Cx*Cy + 2*D(3,z,y,x)*Cx*Cz + 2*D(5,z,y,x)*Cy*Cz
			Dhessianvec = [D(1,z,y,x), D(2,z,y,x), D(3,z,y,x), D(4,z,y,x), D(5,z,y,x), D(6,z,y,x)];
			% Dhmat hessian mat
			DHmat = hessianvaluetomat(Dhessianvec);
			w1mat = Cmat' * DHmat * Cmat - 1/(F(z,y,x)^2);
			% w2 = -2*D(1,z,y,x)*Cx*Kx - 2*D(4,z,y,x)*Cy*Ky - 2*D(6,z,y,x)*Cz*Kz - 2*D(2,z,y,x)*Cx*Ky - 2*D(2,z,y,x)*Cy*Kx - 2*D(3,z,y,x)*Cx*Kz - 2*D(3,z,y,x)*Cz*Kx - 2*D(5,z,y,x)*Cy*Kz - 2*D(5,z,y,x)*Cz*Ky
			w2mat = -2 * Kvec' * DHmat * Cmat;
			% w3 = D(1,z,y,x)*(Kx^2) + D(4,z,y,x)*(Ky^2) + D(6,z,y,x)*(Kz^2) + 2*D(2,z,y,x)*Kx*Ky + 2*D(3,z,y,x)*Kx*Kz + 2*D(5,z,y,x)*Ky*Kz - 1/(F(z,y,x)^2)
			w3mat = Kvec' * DHmat * Kvec;
			[R, flag_imag] = afmfind_roots_indic(100 * w1mat, 100 * w2mat, 100 * w3mat);
			if(~flag_imag && R(1) > 0)
				% dT(1) = p(1)*(R(1) - Ta)/d_a + p(2)*(R(1) - Tb)/d_b + p(3)*(R(1) - Tc)/d_c;
				% dT(2) = q(1)*(R(1) - Ta)/d_a + q(2)*(R(1) - Tb)/d_b + q(3)*(R(1) - Tc)/d_c;
				% dT(3) = r(1)*(R(1) - Ta)/d_a + r(2)*(R(1) - Tb)/d_b + r(3)*(R(1) - Tc)/d_c;
				dTVec = pqrmat *  (R(1) - [Ta; Tb; Tc]).*[1/d_a; 1/d_b; 1/d_c];
				% dT
				% dTVec
				% d(1) = (D(1,z,y,x)*dTVec(1) + D(2,z,y,x)*dTVec(2) + D(3,z,y,x)*dTVec(3));
				% d(2) = (D(2,z,y,x)*dTVec(1) + D(4,z,y,x)*dTVec(2) + D(5,z,y,x)*dTVec(3));
				% d(3) = (D(3,z,y,x)*dTVec(1) + D(5,z,y,x)*dTVec(2) + D(6,z,y,x)*dTVec(3));
				dmat = DHmat * dTVec;
				% dir1 = (dmat(1)*(a(2)*b(3) - a(3)*b(2)) - dmat(2)*(a(1)*b(3) - a(3)*b(1)) + dmat(3)*(a(1)*b(2) - a(2)*b(1)))
				dir1mat = dot(dmat', cross(a, b));  
				% dir2 = (dmat(1)*(b(2)*c(3) - b(3)*c(2)) - dmat(2)*(b(1)*c(3) - b(3)*c(1)) + dmat(3)*(b(1)*c(2) - b(2)*c(1)))
				dir2mat = dot(dmat', cross(b, c));
				% dir3 = (dmat(1)*(c(2)*a(3) - c(3)*a(2)) - dmat(2)*(c(1)*a(3) - c(3)*a(1)) + dmat(3)*(c(1)*a(2) - c(2)*a(1)))
				dir3mat = dot(dmat', cross(c, a));
				conditionone = (afmsign(dir1mat)<0 && afmsign(dir2mat)<=0 && afmsign(dir3mat)<=0);
				conditiontwo = (afmsign(dir1mat)<=0 && afmsign(dir2mat)<=0 && afmsign(dir3mat)<0);
				conditionthree = (afmsign(dir1mat)<=0 && afmsign(dir2mat)<0 && afmsign(dir3mat)<=0);
				conditionfour = true; 
				if(conditionone || conditiontwo || conditionthree && conditionfour)
					temp = afmmin(temp, R(1));
					flag1 = true;
				end
			end
			% R
			if(~flag_imag && R(2) > 0)
				% dT(1) = p(1)*(R(2) - Ta)/d_a + p(2)*(R(2) - Tb)/d_b + p(3)*(R(2) - Tc)/d_c;
				% dT(2) = q(1)*(R(2) - Ta)/d_a + q(2)*(R(2) - Tb)/d_b + q(3)*(R(2) - Tc)/d_c;
				% dT(3) = r(1)*(R(2) - Ta)/d_a + r(2)*(R(2) - Tb)/d_b + r(3)*(R(2) - Tc)/d_c
				dTVec = pqrmat *  (R(2) - [Ta; Tb; Tc]).*[1/d_a; 1/d_b; 1/d_c];
				% dT
				% dTVec
				% d(1) = (D(1,z,y,x)*dTVec(1) + D(2,z,y,x)*dTVec(2) + D(3,z,y,x)*dTVec(3));
				% d(2) = (D(2,z,y,x)*dTVec(1) + D(4,z,y,x)*dTVec(2) + D(5,z,y,x)*dTVec(3));
				% d(3) = (D(3,z,y,x)*dTVec(1) + D(5,z,y,x)*dTVec(2) + D(6,z,y,x)*dTVec(3));
				dmat = DHmat * dTVec;
				% dir1 = (dmat(1)*(a(2)*b(3) - a(3)*b(2)) - dmat(2)*(a(1)*b(3) - a(3)*b(1)) + dmat(3)*(a(1)*b(2) - a(2)*b(1)))
				dir1mat = dot(dmat', cross(a, b));  
				% dir2 = (dmat(1)*(b(2)*c(3) - b(3)*c(2)) - dmat(2)*(b(1)*c(3) - b(3)*c(1)) + dmat(3)*(b(1)*c(2) - b(2)*c(1)))
				dir2mat = dot(dmat', cross(b, c));
				% dir3 = (dmat(1)*(c(2)*a(3) - c(3)*a(2)) - dmat(2)*(c(1)*a(3) - c(3)*a(1)) + dmat(3)*(c(1)*a(2) - c(2)*a(1)))
				dir3mat = dot(dmat', cross(c, a));
				conditionone = (afmsign(dir1mat)<0 && afmsign(dir2mat)<=0 && afmsign(dir3mat)<=0);
				conditiontwo = (afmsign(dir1mat)<=0 && afmsign(dir2mat)<=0 && afmsign(dir3mat)<0);
				conditionthree = (afmsign(dir1mat)<=0 && afmsign(dir2mat)<0 && afmsign(dir3mat)<=0);
				conditionfour = true; 
				if(conditionone || conditiontwo || conditionthree && conditionfour)
					temp = afmmin(temp, R(2));
					flag2 = true;
				end
			end

		end
		
		% if found roots are not useful. both flags should be false.
		% fprintf('the value of flag1 is : %d the value of flag2 is : %d\n', flag1, flag2);
		if( (~flag1 && ~flag2) || 1)
			% all interest points are ok.
			% fprintf('the value of flaga is : %d the value of flagb is : %d the value of flagc is : %d\n', flaga, flagb, flagc);
			% flaga = false; flagb = true; flagc = false; % this line should be deleted in the future !!!!!!!
			if(flaga && flagb && flagc) 
				temp = afmmin(temp, afmminimize_Analytic2(Ta,Tb,x,y,z,D,a,b,F(z,y,x),d_a,d_b));
				temp = afmmin(temp, afmminimize_Analytic2(Ta,Tc,x,y,z,D,a,c,F(z,y,x),d_a,d_c));
				temp = afmmin(temp, afmminimize_Analytic2(Tb,Tc,x,y,z,D,b,c,F(z,y,x),d_b,d_c));
				fprintf('true true true\n');
			elseif(~flaga && flagb && flagc) % only two are ok.
				temp = afmmin(temp, afmminimize_Analytic2(Tb,Tc,x,y,z,D,b,c,F(z,y,x),d_b,d_c));
				fprintf('false true true\n');
			elseif(flaga && ~flagb && flagc) % only two are ok.
				temp = afmmin(temp, afmminimize_Analytic2(Ta,Tc,x,y,z,D,a,c,F(z,y,x),d_a,d_c));
			elseif(flaga && flagb && ~flagc) % only two are ok.
				temp = afmmin(temp, afmminimize_Analytic2(Ta,Tb,x,y,z,D,a,b,F(z,y,x),d_a,d_b));
			elseif(flaga && ~flagb && ~flagc) % only one is ok.
				vga = afmgroup_velocity(x,y,z,D,a',F(z,y,x));
				temp = afmmin(temp, Ta+d_a/vga);
			elseif(~flaga && flagb && ~flagc) % only one is ok.
				vgb = afmgroup_velocity(x,y,z,D,b',F(z,y,x));
				temp = afmmin(temp, Tb+d_b/vgb);
			elseif(~flaga && ~flagb && flagc) % only one is ok.
				vgc = afmgroup_velocity(x,y,z,D,c',F(z,y,x));	
				temp = afmmin(temp, Tc+d_c/vgc);
			end
		end

	end
	% temp = 555 % it will be removed in the future !!!!!!
	if(temp ~= 5000)
		if( Ttag(z,y,x) == 175 ) % point is already in the good trial list. UPDATE.
			viter = value_iter{z,y,x,1};
			firstkey = viter{1};
			if( firstkey > temp )
				Tvalue(z,y,x) = temp;
				index = afmsub2ind(afmSize, x, y, z);
				% fm_map->trial.erase( viter );
				remove(trial, firstkey);
				% value_iter[z][y][x][0] = fm_map->trial.insert( make_pair(temp,index) );
				trial(temp) = index;				
				value_iter{z,y,x,1} = trial;
				% fprintf('Stage One\n');
			end
		elseif( Ttag(z,y,x) == 125 ) % point is not in the good list but is in the bad list. REMOVE AND ADD TO THE GOOD LIST.
			if( Tvalue(z,y,x) >= temp )
				% removing from the bad one.
				viter = value_iter{z,y,x,2};	
				% fm_map->trialC.erase( viter );
				firstkey = viter{1};
				remove(trialC, firstkey);
				% adding to the good trial list
				index = afmsub2ind(afmSize, x,y,z)
				% value_iter[z,y,x,0] = fm_map->trial.insert(make_pair(temp,index)); %  write the address.
				trialC(temp) = index;				
				value_iter{z,y,x,1} = trial;
				Tvalue(z,y,x) = temp;
				Ttag(z,y,x) = 175; % making good trial tag.
				% fprintf('Stage Two\n');
			end
		else % it is not in both lists
			index = afmsub2ind(afmSize, x,y,z)
			% value_iter[z][y][x][0] = fm_map->trial.insert( make_pair(temp,index) ); // write the address.
			% This line will be revised in the future
            temp = real(temp);
            trial(temp) = index;
			value_iter{z,y,x,1} = trial;
			% save('value_iter.mat', 'value_iter');
			Tvalue(z,y,x) = temp;
			Ttag(z,y,x) = 175; % good list trial number
			% fprintf('Stage Three\n');
		end
	end
	% temp2 = 555; %this line will be removed in the future 
	if(temp2 ~= 5000)
		% distance to trial lists.
		if( Ttag(z,y,x) == 175 ) % point is alread in the good trial list. CHECK-REMOVE-UDPATE.
			if(Tvalue(z,y,x) > temp2) % found a better candidate.
				% removing from the good trial list.
				viter = value_iter{z,y,x,1};
				% fm_map->trial.erase(viter);
				firstkey = viter{1};
				remove(trial, firstkey);
				% adding to the bad trial list.
				index = afmsub2ind(afmSize,x,y,z);
				% value_iter[z][y][x][1] = fm_map->trialC.insert( make_pair(temp2,index) ); % write the address
				trialC(temp2) = index;
				value_iter{z,y,x,2} = trialC;				
				Tvalue(z,y,x) = temp2;
				Ttag(z,y,x) = 125; % making bad trial tag.
			end
		elseif( Ttag(z,y,x) == 125 ) % not in the good list but in the bad list.
			viter = value_iter{z,y,x,2};
			firstkey = viter{1};
			if(firstkey > temp2)
				Tvalue(z,y,x) = temp2;
				index = afmsub2ind(afmSize,x,y,z);
				% fm_map->trialC.erase(viter);
				remove(trialC, firstkey);
				% value_iter[z][y][x][1] = fm_map->trialC.insert( make_pair(temp2,index) );
				trialC(temp2) = index;
				value_iter{z,y,x,2} = trialC;	
			end
		else % it is not in both lists
			index = afmsub2ind(afmSize, x,y,z);
			% value_iter[z,y,x,1] = fm_map->trialC.insert( make_pair(temp2,index) ); % write the address
			trialC(temp2) = index;
			value_iter{z,y,x,2} = trialC;	
			Tvalue(z,y,x) = temp2;
			Ttag(z,y,x) = 125; %  bad list trial number
		end
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
	% Tvalue = [];
	% Ttag = [];
	% value_iter = [];
